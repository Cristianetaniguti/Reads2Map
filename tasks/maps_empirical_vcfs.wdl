version 1.0

import "../structs/maps_empirical_vcfsS.wdl"

workflow maps_empirical_vcfs{

input {
    Infos infos
    Vcfs vcfs
}

call BamCounts4Onemap{
    input:
      sampleName = infos.names,
      freebayes_counts = vcfs.freebayes_counts,
      gatk_counts = vcfs.gatk_counts
}

#   Array[String] methods = ["gatk", "freebayes"]
#   Array[File] vcfs = [VcftoolsApplyFilters.gatkVCF_F, VcftoolsApplyFilters.freebayesVCF_F]
#   Array[Pair[String, File]] program_and_vcf = zip(methods, vcfs)

#   scatter (vcf in program_and_vcf){
#     call CreateMaps{
#       input:
#         methodName = vcf.left,
#         vcf_file   = vcf.right,
#         freebayes_ref_depth = BamCounts4Onemap.freebayes_ref_bam,
#         freebayes_alt_depth = BamCounts4Onemap.freebayes_alt_bam,
#         gatk_ref_depth = BamCounts4Onemap.gatk_ref_bam,
#         gatk_alt_depth = BamCounts4Onemap.gatk_alt_bam,
#         gatk_example_alleles = BamCounts4Onemap.gatk_example_alleles,
#         freebayes_example_alleles = BamCounts4Onemap.freebayes_example_alleles,
#         cross = infos.cross
#     }
#   }

#   call CreateTables{
#     input:
#         gatk_ref_depth            = CalculateVcfMetrics.gatk_ref_depth,
#         gatk_ref_depth_bam        = BamCounts4Onemap.gatk_ref_bam,
#         gatk_alt_depth            = CalculateVcfMetrics.gatk_alt_depth,
#         gatk_alt_depth_bam        = BamCounts4Onemap.gatk_alt_bam,
#         freebayes_ref_depth_bam   = BamCounts4Onemap.freebayes_ref_bam,
#         freebayes_alt_depth_bam   = BamCounts4Onemap.freebayes_alt_bam,
#         freebayes_ref_depth       = CalculateVcfMetrics.freebayes_ref_depth,
#         freebayes_alt_depth       = CalculateVcfMetrics.freebayes_alt_depth,
#         map_df                    = CreateMaps.map_df,
#         map_GQ                    = CreateMaps.map_GQ,
#         map_polyrad               = CreateMaps.map_polyrad,
#         map_supermassa            = CreateMaps.map_supermassa,
#         map_updog                 = CreateMaps.map_updog,
#         map_bam_polyrad           = CreateMaps.map_bam_polyrad,
#         map_bam_supermassa        = CreateMaps.map_bam_supermassa,
#         map_bam_updog             = CreateMaps.map_bam_updog,
#         error_info_df             = CreateMaps.error_info_df,
#         error_info_GQ             = CreateMaps.error_info_GQ,
#         error_info_updog          = CreateMaps.error_info_updog,
#         error_info_polyrad        = CreateMaps.error_info_polyrad,
#         error_info_supermassa     = CreateMaps.error_info_supermassa,
#         error_info_bam_updog      = CreateMaps.error_info_bam_updog,
#         error_info_bam_polyrad    = CreateMaps.error_info_bam_polyrad,
#         error_info_bam_supermassa = CreateMaps.error_info_bam_supermassa,
#         filters_dfAndGQ           = CreateMaps.filters_dfAndGQ,
#         filters_polyrad           = CreateMaps.filters_polyrad,
#         filters_supermassa        = CreateMaps.filters_supermassa,
#         filters_updog             = CreateMaps.filters_updog,
#         filters_bam_polyrad       = CreateMaps.filters_bam_polyrad,
#         filters_bam_supermassa    = CreateMaps.filters_bam_supermassa,
#         filters_bam_updog         = CreateMaps.filters_bam_updog
#     }

#   output {
#     File data1 = CreateTables.data1
#     File data2 = CreateTables.data2
#     File data3 = CreateTables.data3
#     File data4 = CreateTables.data4
#     File data5 = CalculateVcfMetrics.data5
#   }
}


# This task convert the output from BamCounts to the depths input for onemap
task BamCounts4Onemap{
  input{
    Array[File] freebayes_counts
    Array[File] gatk_counts
    Array[String] sampleName
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(R.utils)
      system("cp ~{sep=" "  freebayes_counts} .")
      system("cp ~{sep=" "  gatk_counts} .")
      names <- c("~{sep=" , "  sampleName}")
      names <- unlist(strsplit(names, split = " , "))

      methods <- c("gatk", "freebayes")

      for(method in methods){

          file.counts <- read.table(paste0(names[1],"_", method,"_counts.tsv"), skip = 3, header=T, stringsAsFactors = F)

          ref_depth_matrix2 <- alt_depth_matrix2  <- matrix(NA, nrow = dim(file.counts)[1], ncol = length(names))

          for(j in 1:length(names)){
            ## From picard tool

            file.counts <- read.table(paste0(names[j],"_", method,"_counts.tsv"), skip = 3, header=T, stringsAsFactors = F)

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
          colnames(ref_depth_matrix2) <- colnames(alt_depth_matrix2) <- names

          alleles <- data.frame(file.counts[,1],file.counts[,2], ref_allele, alt_allele)
          write.table(alleles, file = paste0(method,"_example4ref_alt_alleles.txt"), col.names = F, row.names = F)

          write.table(ref_depth_matrix2, file = paste0(method,"_ref_depth_bam.txt"), quote=F, row.names=T, sep="\t", col.names=T)
          write.table(alt_depth_matrix2, file = paste0(method,"_alt_depth_bam.txt"), quote=F, row.names=T, sep="\t", col.names=T)
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