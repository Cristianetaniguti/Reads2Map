version 1.0

task ProcessRadTags {
    input {
      String enzyme
      String? enzyme2
      Array[File] fq_files
      File? barcodes
    }

    Int disk_size = ceil(size(fq_files, "GiB") * 2)
    Int memory_size = 6000

    command <<<
      mkdir raw process_radtags_results
      mv ~{sep=" " fq_files} raw

      process_radtags -p raw/ -o process_radtags_results/ \
                      ~{"-b " + barcodes} \
                      --renz_1 ~{enzyme} ~{"--renz_2 " + enzyme2} \
                      -r -c -q -w 0.5

    >>>

    runtime {
      docker:"cristaniguti/stacks:0.0.1"
      cpu:1
      # Cloud
      memory:"~{memory_size} MiB"
      disks:"local-disk " + disk_size + " HDD"
      # Slurm
      job_name: "ProcessRadTags"
      mem:"~{memory_size}M"
      time:"10:00:00"
    }

    output {
      Array[File] seq_results = glob("process_radtags_results/*.fq.gz")
    }
}

task CreatePopMapFile {
   input {
      Array[File] bams
   }
   
   command <<<
    R --vanilla --no-save <<RSCRIPT

      file <- "~{sep=',' bams}"
      files <- unlist(strsplit(file, ","))
      file_names <- gsub(".bam","",basename(files))

      pop_file <- data.frame(file_names, "pop1")
      write.table(pop_file, file = "pop_map.tsv", 
                  sep = "\t", quote = FALSE, row.names = FALSE, 
                  col.names = FALSE)

    RSCRIPT
   >>>

   runtime {
      docker:"cristaniguti/r-samtools:latest"
      cpu:1
      # Cloud
      memory:"100 MiB"
      disks:"local-disk 1 GiB HDD"
      # Slurm
      job_name: "CreatePopMapFile"
      mem:"100M"
      time:"10:00:00"
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Reads a SAM file to determine the potential positions of Tags against the reference genome. Identifies SNPs from the aligned tags"
    }

    output {
      File pop_map = "pop_map.tsv"
    }
}

task RefMap {
  input {
    Array[File] bams
    File pop_map
    Int max_cores
  }

    Int disk_size = ceil(size(bams, "GiB") * 2)
    Int memory_size = 6000

  command <<<

    mkdir aligned stacks
    ln -s ~{sep=" " bams} aligned

    gstacks -I aligned/ -M ~{pop_map} -O stacks/ -t ~{max_cores}
        
    populations -P stacks/ -M ~{pop_map} --vcf -t ~{max_cores}

  >>>

   runtime {
      docker:"cristaniguti/stacks:0.0.1"
      cpu:1
      # Cloud
      memory:"~{memory_size} MiB"
      disks:"local-disk " + disk_size + " HDD"
      # Slurm
      job_name: "ProcessRadTags"
      mem:"~{memory_size}M"
      time:"10:00:00"
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Reads a SAM file to determine the potential positions of Tags against the reference genome. Identifies SNPs from the aligned tags"
    }

    output {
      File stacks_vcf = "stacks/populations.snps.vcf"
      File stacks_multiallelics = "stacks/populations.haps.vcf"
    }
}
