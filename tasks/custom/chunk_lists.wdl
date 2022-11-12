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

