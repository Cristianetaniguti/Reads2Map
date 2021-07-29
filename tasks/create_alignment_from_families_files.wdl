version 1.0

import "alignment.wdl" as alg

workflow CreateAlignmentFromFamilies {
    input {
        File families_info
        Reference references
        Int max_cores
        String rm_dupli
        Int chunk_size
    }

    call SepareChunks {
        input:
            families_info=families_info,
            chunk_size = chunk_size
    }

    scatter (chunk in SepareChunks.chunks) {

        Array[Array[String]] sample_file = read_tsv(chunk)

        call alg.RunBwaAlignment {
            input:
                sampleName  = sample_file[1],
                reads       = sample_file[0],
                libraries   = sample_file[2],
                references  = references,
                max_cores   = max_cores,
                rm_dupli    = rm_dupli
        }
    }

    output {
        Array[File] bam = flatten(RunBwaAlignment.bam)
        Array[File] bai = flatten(RunBwaAlignment.bai)
        Array[File] dup_metrics = flatten(RunBwaAlignment.dup_metrics)
    }
}

task SepareChunks {
    input {
        File families_info
        Int chunk_size
    }

    command <<<
        R --vanilla --no-save <<RSCRIPT
            library(rjson)

            df <- read.table("~{families_info}")
            split_df <- split.data.frame(df, df[,2])

            n_chunk <- as.integer(length(split_df)/~{chunk_size})
            chunk_temp <- rep(1:n_chunk, each=~{chunk_size})
            chunk <- c(chunk_temp, rep(n_chunk+1, length(split_df) - length(chunk_temp)))
            chunk_json <- split(split_df, chunk)

            for(i in 1:length(chunk_json)){
                df <- do.call(rbind, unlist(chunk_json[i], recursive = F))
                df <- t(df)
                write.table(df, file = paste0("chunk_",i, ".txt"), quote = F, col.names = F, row.names = F, sep="\t")
            }

        RSCRIPT

    >>>

    runtime {
        job_name: "SepareChunksIndividuals"
        docker: "cristaniguti/r-samtools"
        node:"--nodes=1"
        mem:"--mem=1G"
        tasks:"--ntasks=1"
        time:"00:05:00"
    }

    output {
        Array[File] chunks = glob("chunk*")
    }
}