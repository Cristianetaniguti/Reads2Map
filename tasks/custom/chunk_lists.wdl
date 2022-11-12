version 1.0

task SepareChunks {
    input {
        File families_info
        Int chunk_size
    }

    Int disk_size = ceil(size(families_info, "GiB") * 2)
    Int memory_size = 1000

    command <<<
        R --vanilla --no-save <<RSCRIPT
            df <- read.table("~{families_info}")
            split_df <- split.data.frame(df, df[,2])

            n_chunk <- as.integer(length(split_df)/~{chunk_size})
            chunk_temp <- rep(1:n_chunk, each=~{chunk_size})
            chunk <- c(chunk_temp, rep(n_chunk+1, length(split_df) - length(chunk_temp)))
            chunk_sep <- split(split_df, chunk)

            for(i in 1:length(chunk_sep)){
                df <- do.call(rbind, unlist(chunk_sep[i], recursive = F))
                df <- t(df)
                write.table(df, file = paste0("chunk_",i, ".txt"), quote = F, col.names = F, row.names = F, sep="\t")
            }

        RSCRIPT

    >>>

    runtime {
        docker: "cristaniguti/reads2map:0.0.1"
        cpu:1
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "SepareChunksIndividuals"
        mem:"~{memory_size}M"
        time:"00:10:00"
    }

    meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Split the fastq files into chunks to be aligned in parallel in the next task."
    }

    output {
        Array[File] chunks = glob("chunk*")
    }
}


task SepareChunksInputList {
  input {
      Array[File] fastqs
      Int chunk_size
  }

  Int disk_size = ceil(size(fastqs, "GiB") * 2)
  Int memory_size = 1000

  command <<<
        R --vanilla --no-save <<RSCRIPT

            files <- c("~{sep="," fastqs}")
            files <- unlist(strsplit(files, split = ","))

            n_chunk <- floor(length(files)/~{chunk_size})

            chunk_temp <- rep(1:n_chunk, each = ~{chunk_size})
            chunk <- c(chunk_temp, rep(n_chunk+1, length(files) - length(chunk_temp)))

            chunk_sep <- split(files, chunk)

            for(i in 1:length(chunk_sep)){
              write.table(chunk_sep[[i]], file = paste0("chunk_",i, ".txt"), quote = F, col.names = F, row.names = F, sep="\t")
            }

        RSCRIPT

  >>>

  runtime {
      docker: "cristaniguti/reads2map:0.0.1"
      cpu:1
      # Cloud
      memory:"~{memory_size} MiB"
      disks:"local-disk " + disk_size + " HDD"
      # Slurm
      job_name: "SepareChunksIndividuals"
      mem:"~{memory_size}M"
      time:"00:10:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Split the simulated fastq files into chunks to be aligned in parallel in the next task."
  }

  output {
    Array[File] chunks = glob("chunk*")
  }
}


task CreateChunks {
  input {
    Array[String] bams
    Array[String] bams_index
    File reference_fasta
    Int chunk_size
  }

  Int disk_size = ceil(size(reference_fasta, "GiB") + 2)

  command <<<
    set -e
    for i in ~{sep=" " bams}; do echo $i >> lof_bams.txt; done
    for i in ~{sep=" " bams_index}; do echo $i >> lof_bais.txt; done

    split -l ~{chunk_size} lof_bams.txt bams.
    split -l ~{chunk_size} lof_bais.txt bais.

    cat ~{reference_fasta} | grep '>' | tr '\n' ',' | sed '$ s/.$//' | sed 's/,/ \n/g' | sed 's/>//g' > intervals.txt
  >>>

  runtime {
    docker: "ubuntu:20.04"
    cpu: 1
    # Cloud
    memory:"1000 MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "CreateChunks"
    mem:"1G"
    time:"00:05:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Split the the samples BAM alignment files in chunks."
  }

  output {
    Array[File] bams_chunks = glob("bams.*")
    Array[File] bais_chunks = glob("bais.*")
    File interval_list = "intervals.txt"
  }
}


task CreateChunksFromGatkPolyploid {  # TODO: We can use SepareChunks and remove this one that came from gatk_polyploid.wdl
  input {
    Array[String] bams
    Array[String] bams_index
    File reference_fasta
    Int chunk_size
  }

  Int disk_size = ceil(size(reference_fasta, "GiB") + 2)

  command <<<
    set -e
    for i in ~{sep=" " bams}; do echo $i >> lof_bams.txt; done
    for i in ~{sep=" " bams_index}; do echo $i >> lof_bais.txt; done

    split -l ~{chunk_size} lof_bams.txt bams.
    split -l ~{chunk_size} lof_bais.txt bais.

    cat ~{reference_fasta} | grep '>' | tr '\n' ',' | sed '$ s/.$//' | sed 's/,/ \n/g' | sed 's/>//g' > intervals.txt
  >>>

  runtime {
    docker: "ubuntu:20.04"
    cpu: 1
    # Cloud
    memory:"1000 MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "CreateChunks"
    mem:"1G"
    time:"00:05:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Split the the samples BAM alignment files in chunks."
  }

  output {
    Array[File] bams_chunks = glob("bams.*")
    Array[File] bais_chunks = glob("bais.*")
    File interval_list = "intervals.txt"
  }
}


task SepareChunksFromGatkPolyploid {  # TODO: We can use SepareChunks and remove this one that came from gatk_polyploid.wdl
    input {
        File families_info
        Int chunk_size
    }

    command <<<
        R --vanilla --no-save <<RSCRIPT
            df <- read.table("~{families_info}")
            split_df <- split.data.frame(df, df[,2])

            n_chunk <- as.integer(length(split_df)/~{chunk_size})
            chunk_temp <- rep(1:n_chunk, each=~{chunk_size})
            chunk <- c(chunk_temp, rep(n_chunk+1, length(split_df) - length(chunk_temp)))
            chunk_sep <- split(split_df, chunk)

            for(i in 1:length(chunk_sep)){
                df <- do.call(rbind, unlist(chunk_sep[i], recursive = F))
                df <- t(df)
                write.table(df, file = paste0("chunk_",i, ".txt"), quote = F, col.names = F, row.names = F, sep="\t")
            }

        RSCRIPT

    >>>

    runtime {
        docker: "cristaniguti/r-samtools:latest"
        cpu: 1
        # Cloud
        memory:"200 MiB"
        disks:"1 HDD"
        # Slurm
        job_name: "SepareChunksIndividuals"
        mem:"1G"
        time:"00:05:00"
    }

    meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Split the fastq files into chunks to be aligned in parallel in the next task."
    }

    output {
        Array[File] chunks = glob("chunk*")
    }
}


task SepareChunksBed {
    input {
        File bed_file
        Int n_nodes
    }

    Int disk_size = ceil(size(bed_file, "GiB") * 1.5)
    Int memory_size = 1000

    command <<<
        R --vanilla --no-save <<RSCRIPT
            library(vroom)
            bed <- vroom("~{bed_file}")

            chunk_size <- floor(dim(bed)[1]/~{n_nodes})

            chunk_temp <- rep(1:~{n_nodes}, each = chunk_size)
            chunk <- c(chunk_temp, rep(~{n_nodes}+1, dim(bed)[1] - length(chunk_temp)))

            chunk_sep <- split.data.frame(bed, chunk)

            for(i in 1:length(chunk_sep)){
                write.table(chunk_sep[[i]], file = paste0("chunk_",i, ".txt"), quote = F, col.names = F, row.names = F, sep="\t")
            }
        RSCRIPT
    >>>

    runtime {
        docker: "cristaniguti/reads2map:0.0.1"
        cpu: 1
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "SepareChunksBed"
        mem:"~{memory_size}M"
        time:"00:05:00"
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Split BED file rows in batches according to defined number of nodes."
    }

    output {
        Array[File] chunks = glob("chunk*")
    }
}
