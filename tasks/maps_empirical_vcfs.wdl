version 1.0

import "../structs/maps_empirical_vcfsS.wdl"

workflow maps_empirical_vcfs{

  input {
      Infos infos
      Vcfs vcfs
  }
  
  Array[File] tsv_files = read_lines(infos.tsv_counts)

  call BamCounts4OneMap{
      input:
        tsv_counts = tsv.files
  }
}

task BamCounts4OneMap{
  input{
    File tsv_counts
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(R.utils)
      tsv_files <- read.table("~{tsv_counts}", stringsAsFactors = F)
      tsv_files <- tsv_files[,1]
      freebayes_counts <- tsv_files[grep("freebayes", tsv_files)]
      gatk_counts <- tsv_files[grep("gatk", tsv_files)]
      counts <- list(freebayes_counts, gatk_counts)

      methods <- c("freebayes", "gatk")

      # BamCounts4OneMap
      for(method in 1:2){

        file.counts <- read.table(counts[[method]][1], skip = 1448, header=T, stringsAsFactors = F)

        ref_depth_matrix2 <- alt_depth_matrix2  <- matrix(NA, nrow = dim(file.counts)[1], ncol = length(counts[[method]]))

        for(j in 1:length(counts[[method]])){
          ## From picard tool

          file.counts <- read.table(counts[[method]][j], skip = 1448, header=T, stringsAsFactors = F)

          ref_depth_matrix2[,j] <- file.counts[,3]
          alt_depth_matrix2[,j] <- file.counts[,4]

          if(j == 1){
            ref_allele <- file.counts[,5]
            alt_allele <- file.counts[,6]
          } else {
            idx.ref <- which(ref_allele == "N")
            idx.alt <- which(alt_allele == "N")
            if(length(idx.ref)!=0){
              ref_allele[idx.ref] <- file.counts[idx.ref,5]
            }
            if(length(idx.alt)!=0){
              alt_allele[idx.alt] <- file.counts[idx.alt,6]
            }
          }

        }

        rownames(ref_depth_matrix2) <- rownames(alt_depth_matrix2) <- paste0(file.counts[,1],"_", file.counts[,2])
        
        alleles <- data.frame(file.counts[,1],file.counts[,2], ref_allele, alt_allele)
        write.table(alleles, file = paste0(methods[method],"_example4ref_alt_alleles.txt"), col.names = F, row.names = F)

        write.table(ref_depth_matrix2, file = paste0(methods[method],"_ref_depth_bam.txt"), quote=F, row.names=T, sep="\t", col.names=T)
        write.table(alt_depth_matrix2, file = paste0(methods[method],"_alt_depth_bam.txt"), quote=F, row.names=T, sep="\t", col.names=T)
      }

    RSCRIPT

  >>>

  runtime{
    docker:"taniguti/onemap"
  }

  output{
    File gatk_ref_bam      = "gatk_ref_depth_bam.txt"
    File gatk_alt_bam      = "gatk_alt_depth_bam.txt"
    File freebayes_ref_bam = "freebayes_ref_depth_bam.txt"
    File freebayes_alt_bam = "freebayes_alt_depth_bam.txt"
    File gatk_example_alleles    = "gatk_example4ref_alt_alleles.txt"
    File freebayes_example_alleles    = "freebayes_example4ref_alt_alleles.txt"

  }
}