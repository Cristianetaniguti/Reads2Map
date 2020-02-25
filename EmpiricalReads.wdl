version 1.0

import "structs/empiricalS.wdl"
import "tasks/create_alignment_from_families_files.wdl" as fam
import "tasks/gatk_genotyping.wdl" as gatk
import "tasks/freebayes_genotyping.wdl" as freebayes
import "tasks/utils.wdl" as utils


workflow EmpiricalReads {

  input {
    Dataset dataset
    ReferenceFasta references
  }

  call fam.CreateAlignmentFromFamilies {
    input:
      families_info=dataset.samples_info,
      references=references
  }

  call gatk.GatkGenotyping {
    input:
      alignments=CreateAlignmentFromFamilies.alignments,
      references=references,
      program="gatk"
  }

  call freebayes.FreebayesGenotyping {
    input:
      alignments=CreateAlignmentFromFamilies.alignments,
      bam=CreateAlignmentFromFamilies.bam,
      bai=CreateAlignmentFromFamilies.bai,
      references=references,
      program="freebayes"
  }

  call utils.BamCounts4Onemap {
    input:
      sampleName=CreateAlignmentFromFamilies.names,
      freebayes_counts=FreebayesGenotyping.counts,
      gatk_counts=GatkGenotyping.counts
  }

  Array[String] methods = ["gatk", "freebayes"]
  Array[File] vcfs = [GatkGenotyping.vcf, FreebayesGenotyping.vcf]
  Array[Pair[String, File]] program_and_vcf = zip(methods, vcfs)

  scatter (vcf in program_and_vcf){
    call avalVCFs {
      input:
      methodName=vcf.left,
      vcf_file=vcf.right,
      freebayes_ref_depth=BamCounts4Onemap.freebayes_ref_bam,
      freebayes_alt_depth=BamCounts4Onemap.freebayes_alt_bam,
      gatk_ref_depth=BamCounts4Onemap.gatk_ref_bam,
      gatk_alt_depth=BamCounts4Onemap.gatk_alt_bam,
      parent1=dataset.parent1,
      parent2=dataset.parent2,
      gatk_example_alleles=BamCounts4Onemap.gatk_example_alleles,
      freebayes_example_alleles=BamCounts4Onemap.freebayes_example_alleles
    }
  }

  call JointDatas{
    input:
      depths_GQ_vcf_rds = avalVCFs.depths_GQ_vcf_rds,
      depths_df_vcf_rds = avalVCFs.depths_df_vcf_rds,
      depths_updog_vcf_rds = avalVCFs.depths_updog_vcf_rds,
      depths_supermassa_vcf_rds = avalVCFs.depths_supermassa_vcf_rds,
      depths_polyrad_vcf_rds = avalVCFs.depths_polyrad_vcf_rds,
      depths_updog_bam_rds = avalVCFs.depths_updog_bam_rds,
      depths_supermassa_bam_rds = avalVCFs.depths_supermassa_bam_rds,
      depths_polyrad_bam_rds = avalVCFs.depths_polyrad_bam_rds,
      times_vcf_df = avalVCFs.times_vcf_df,
      times_vcf_GQ = avalVCFs.times_vcf_GQ,
      times_vcf_updog = avalVCFs.times_vcf_updog,
      times_vcf_supermassa = avalVCFs.times_vcf_supermassa,
      times_vcf_polyrad = avalVCFs.times_vcf_polyrad,
      times_vcf_gusmap = avalVCFs.times_vcf_gusmap,
      times_bam_gusmap = avalVCFs.times_bam_gusmap,
      times_bam_updog = avalVCFs.times_bam_updog,
      times_bam_supermassa = avalVCFs.times_bam_supermassa,
      times_bam_polyrad = avalVCFs.times_bam_polyrad,
      filters_vcf_dfAndGQ = avalVCFs.filters_vcf_dfAndGQ,
      filters_vcf_polyrad = avalVCFs.filters_vcf_polyrad,
      filters_vcf_supermassa = avalVCFs.filters_vcf_supermassa,
      filters_vcf_updog = avalVCFs.filters_vcf_updog,
      filters_bam_polyrad = avalVCFs.filters_bam_polyrad,
      filters_bam_supermassa = avalVCFs.filters_bam_supermassa,
      filters_bam_updog = avalVCFs.filters_bam_updog,
      RData_vcf_df = avalVCFs.RData_vcf_df,
      RData_vcf_GQ = avalVCFs.RData_vcf_GQ,
      RData_vcf_polyrad = avalVCFs.RData_vcf_polyrad,
      RData_vcf_supermassa = avalVCFs.RData_vcf_supermassa,
      RData_vcf_updog = avalVCFs.RData_vcf_updog,
      RData_bam_polyrad = avalVCFs.RData_bam_polyrad,
      RData_bam_supermassa = avalVCFs.RData_bam_supermassa,
      RData_bam_updog = avalVCFs.RData_bam_updog,
      RData_gusmap = avalVCFs.RData_gusmap,
      RData_bam_gusmap = avalVCFs.RData_bam_gusmap,
      map_vcf_df = avalVCFs.map_vcf_df,
      map_vcf_GQ = avalVCFs.map_vcf_GQ,
      map_vcf_polyrad = avalVCFs.map_vcf_polyrad,
      map_vcf_supermassa = avalVCFs.map_vcf_supermassa,
      map_vcf_updog = avalVCFs.map_vcf_updog,
      map_bam_polyrad = avalVCFs.map_bam_polyrad,
      map_bam_supermassa = avalVCFs.map_bam_supermassa,
      map_bam_updog = avalVCFs.map_bam_updog,
      map_gusmap = avalVCFs.map_gusmap,
      map_bam_gusmap = avalVCFs.map_bam_gusmap
    }

  output {
    File data_depths = JointDatas.data_depths
    File data_filters = JointDatas.data_filters
    File data_map = JointDatas.data_map
    File data_times = JointDatas.data_times
    File data_RDatas = JointDatas.data_RDatas
  }
}

task avalVCFs {
  input {
    String methodName
    File vcf_file
    File freebayes_ref_depth
    File freebayes_alt_depth
    File gatk_ref_depth
    File gatk_alt_depth
    String parent1
    String parent2
    File gatk_example_alleles
    File freebayes_example_alleles
  }

 command <<<

   R --vanilla --no-save <<RSCRIPT

   system("cp ~{freebayes_ref_depth} .")
   system("cp ~{freebayes_alt_depth} .")
   system("cp ~{gatk_ref_depth} .")
   system("cp ~{gatk_alt_depth} .")
   system("cp ~{gatk_example_alleles} .")
   system("cp ~{freebayes_example_alleles} .")

   source("/opt/scripts/functions_empirical.R")

   vcf_file <- "~{vcf_file}"
   SNPCall <- "~{methodName}"
   cross <- "outcross"
   max.cores <- 3
   parent1 <- "~{parent1}"
   parent2 <- "~{parent2}"

   # PACKAGES
   library(supermassa4onemap)
   library(onemap)
   library(updog)
   library(reshape2)
   library(vcfR)
   library(doParallel)
   library(GUSMap)
   library(ggplot2)

   ## READING FINAL VCF FROM PIPELINE
   vcf <- read.vcfR(vcf_file)
   df <- onemap_read_vcfR(vcfR.object=vcf,
                         cross= cross,
                         parent1= parent1,
                         parent2= parent2)

   # check depths
   p <- create_depths_profile(onemap.obj = df, vcfR.object = vcf, parent1 = parent1, parent2 = parent2, vcf.par = "AD",
   recovering = FALSE, GTfrom = "vcf", alpha=0.1,
   rds.file = paste0(SNPCall,"_", GenoCall= "df","_",CountsFrom="vcf","_depths.rds"))
   #ggsave(filename = paste0(SNPCall,"_", GenoCall= "df","_",CountsFrom="vcf","_vcf_depths.png"), p)

   ## FILTERS REPORT
   out_name <- paste0(SNPCall, "_filters_vcf_dfAndGQ.txt")
   filters_tab <- create_filters_report(onemap_obj = df, CountsFrom = "vcf", SNPCall=SNPCall, GenoCall="df")
   write_report(filters_tab[[1]], out_name)

   # MAPS REPORT - DF
   create_map_report(input.seq = filters_tab[[2]], CountsFrom = "vcf", SNPCall = SNPCall, GenoCall="df")

   # MAPS REPORT - GQ
   aval.gq <- extract_depth(vcfR.object=vcf,
                           onemap.object=df,
                           vcf.par="GQ",
                           parent1=parent1,
                           parent2=parent2,
                           recovering=FALSE)

   aval.gq <- create_probs(df, genotypes_errors=aval.gq)
   filters_tab <- create_filters_report(aval.gq, CountsFrom = "vcf", SNPCall=SNPCall, GenoCall="df")

   # check depths and errors
   p <- create_depths_profile(onemap.obj = aval.gq, vcfR.object = vcf, parent1 = parent1, parent2 = parent2, vcf.par = "AD",
   recovering = FALSE, GTfrom = "vcf", alpha=0.1,
   rds.file = paste0(SNPCall,"_", GenoCall= "GQ","_",CountsFrom="vcf","_depths.rds"))

   create_map_report(filters_tab[[2]], CountsFrom = "vcf", SNPCall = SNPCall, GenoCall="GQ")

   # OTHER TOOLS
   ## With depths from vcf

   updog.aval <- updog_error(
     vcfR.object=vcf,
     onemap.object=df,
     vcf.par="AD",
     parent1=parent1,
     parent2=parent2,
     recovering=TRUE,
     mean_phred=20,
     cores=max.cores,
     depths=NULL)

   supermassa.aval <- supermassa4onemap::supermassa_error(
     vcfR.object=vcf,
     onemap.object = df,
     vcf.par = "AD",
     parent1 = parent1,
     parent2 = parent2,
     recovering = TRUE,
     mean_phred = 20,
     cores = max.cores,
     depths = NULL)

   polyrad.aval <- polyRAD_error(
     vcf=vcf_file,
     onemap.obj=df,
     parent1=parent1,
     parent2=parent2,
     crosstype=cross)

   metodologies <- list(updog = updog.aval, supermassa = supermassa.aval, polyrad = polyrad.aval)
   for (metodology in names(metodologies)){
     cat(metodology, "\n")
     error_aval <- metodologies[[metodology]]
     ## Filters
     out_name <- paste0(SNPCall, "_filters_vcf_", metodology, ".txt")
     filters_tab <- create_filters_report(error_aval, CountsFrom = "vcf", SNPCall = SNPCall, GenoCall = metodology)
     write_report(filters_tab[[1]], out_name)

     # check depths
     p <- create_depths_profile(onemap.obj = error_aval, vcfR.object = vcf, parent1 = parent1, parent2 = parent2, vcf.par = "AD",
                               recovering = FALSE, GTfrom = "onemap", alpha=0.1,
                               rds.file = paste0(SNPCall,"_", GenoCall= metodology, "_",CountsFrom = "vcf","_depths.rds"))

     #ggsave(filename = paste0(SNPCall,"_", GenoCall= metodology, "_",CountsFrom = "vcf","_onemap_depths.png"), p)

     ## Maps
     create_map_report(input.seq = filters_tab[[2]], CountsFrom = "vcf", SNPCall = SNPCall, GenoCall= metodology)
   }

   ## Depths from bam
   depths.alt <- read.table(paste0(SNPCall, "_alt_depth_bam.txt"), header = T)
   depths.ref <- read.table(paste0(SNPCall, "_ref_depth_bam.txt"), header = T)

   depths <- list("ref"=depths.ref, "alt"=depths.alt)

   updog.aval.bam <- updog_error(
     vcfR.object=vcf,
     onemap.object=df,
     vcf.par="AD",
     parent1=parent1,
     parent2=parent2,
     recovering=TRUE,
     mean_phred=20,
     cores=max.cores,
     depths=depths)

   supermassa.aval.bam <- supermassa_error(
     vcfR.object=vcf,
     onemap.object = df,
     vcf.par = "AD",
     parent1 = parent1,
     parent2 = parent2,
     recovering = TRUE,
     mean_phred = 20,
     cores = max.cores,
     depths = depths)

   if(tail(strsplit(vcf_file, "[.]")[[1]],1) =="gz") {
       vcf.temp <- paste0(SNPCall,".", sample(1000,1), ".vcf")
       system(paste0("zcat ", vcf_file, " > ", vcf.temp))
       vcf_file <- vcf.temp
   }

   new.vcf <- make_vcf(vcf_file, depths, SNPCall)
   new.vcfR <- read.vcfR(new.vcf)

   polyrad.aval.bam <- polyRAD_error(
     vcf=new.vcf,
     onemap.obj=df,
     parent1=parent1,
     parent2=parent2,
     crosstype=cross)

   metodologies <- list(updog = updog.aval.bam, supermassa= supermassa.aval.bam, polyrad=polyrad.aval.bam)
   for (metodology in names(metodologies)){
     cat(metodology, "\n")
     error_aval <- metodologies[[metodology]]
     ## Filters
     out_name <- paste0(SNPCall, "_filters_bam_", metodology, ".txt")
     filters_tab <- create_filters_report(error_aval, CountsFrom = "bam", SNPCall = SNPCall, GenoCall = metodology)
     write_report(filters_tab[[1]], out_name)

     # check depths
     p <- create_depths_profile(onemap.obj = error_aval, vcfR.object = new.vcfR, parent1 = parent1, parent2 = parent2, vcf.par = "AD",
                               recovering = FALSE, GTfrom = "onemap", alpha=0.1,
                               rds.file = paste0(SNPCall,"_", GenoCall= metodology, "_",CountsFrom = "bam","_depths.rds"))
     #ggsave(filename = paste0(SNPCall,"_", GenoCall= metodology, "_",CountsFrom = "bam","_onemap_depths.png"), p)

     ## Maps
     create_map_report(input.seq = filters_tab[[2]], CountsFrom = "bam", SNPCall = SNPCall, GenoCall= metodology)
   }

   ## Gusmap maps
   out_name <- paste0(SNPCall, "_vcf_gusmap_map.txt")
   map_gus <- create_gusmap_report(vcf_file,"gusmap", "vcf", SNPCall, parent1, parent2)
   write_report(map_gus, out_name)

   out_name <- paste0(SNPCall, "_bam_gusmap_map.txt")
   map_gus <- create_gusmap_report(new.vcf, "gusmap", "bam", SNPCall, parent1, parent2)
   write_report(map_gus, out_name)

   RSCRIPT
 >>>

 runtime{
   docker:"taniguti/onemap"
 }

 output{
   File depths_GQ_vcf_rds = "~{methodName}_GQ_vcf_depths.rds"
   File depths_df_vcf_rds = "~{methodName}_df_vcf_depths.rds"
   File depths_updog_vcf_rds = "~{methodName}_updog_vcf_depths.rds"
   File depths_supermassa_vcf_rds = "~{methodName}_supermassa_vcf_depths.rds"
   File depths_polyrad_vcf_rds = "~{methodName}_polyrad_vcf_depths.rds"
   File depths_updog_bam_rds = "~{methodName}_updog_bam_depths.rds"
   File depths_supermassa_bam_rds = "~{methodName}_supermassa_bam_depths.rds"
   File depths_polyrad_bam_rds = "~{methodName}_polyrad_bam_depths.rds"
   File times_vcf_df = "~{methodName}_vcf_df_times.txt"
   File times_vcf_GQ = "~{methodName}_vcf_GQ_times.txt"
   File times_vcf_updog = "~{methodName}_vcf_updog_times.txt"
   File times_vcf_supermassa = "~{methodName}_vcf_supermassa_times.txt"
   File times_vcf_polyrad = "~{methodName}_vcf_polyrad_times.txt"
   File times_vcf_gusmap = "~{methodName}_vcf_gusmap_times.txt"
   File times_bam_gusmap = "~{methodName}_bam_gusmap_times.txt"
   File times_bam_updog = "~{methodName}_bam_updog_times.txt"
   File times_bam_supermassa = "~{methodName}_bam_supermassa_times.txt"
   File times_bam_polyrad = "~{methodName}_bam_polyrad_times.txt"
   File filters_vcf_dfAndGQ = "~{methodName}_filters_vcf_dfAndGQ.txt"
   File filters_vcf_polyrad = "~{methodName}_filters_vcf_polyrad.txt"
   File filters_vcf_supermassa = "~{methodName}_filters_vcf_supermassa.txt"
   File filters_vcf_updog = "~{methodName}_filters_vcf_updog.txt"
   File filters_bam_polyrad = "~{methodName}_filters_bam_polyrad.txt"
   File filters_bam_supermassa = "~{methodName}_filters_bam_supermassa.txt"
   File filters_bam_updog = "~{methodName}_filters_bam_updog.txt"
   File RData_vcf_df = "~{methodName}_vcf_df.RData"
   File RData_vcf_GQ = "~{methodName}_vcf_GQ.RData"
   File RData_vcf_polyrad = "~{methodName}_vcf_polyrad.RData"
   File RData_vcf_supermassa = "~{methodName}_vcf_supermassa.RData"
   File RData_vcf_updog = "~{methodName}_vcf_updog.RData"
   File RData_bam_polyrad = "~{methodName}_bam_polyrad.RData"
   File RData_bam_supermassa = "~{methodName}_bam_supermassa.RData"
   File RData_bam_updog = "~{methodName}_bam_updog.RData"
   File RData_gusmap = "~{methodName}_vcf_gusmap.RData"
   File RData_bam_gusmap = "~{methodName}_bam_gusmap.RData"
   File map_vcf_df = "~{methodName}_vcf_df_map.txt"
   File map_vcf_GQ = "~{methodName}_vcf_GQ_map.txt"
   File map_vcf_polyrad = "~{methodName}_vcf_polyrad_map.txt"
   File map_vcf_supermassa = "~{methodName}_vcf_supermassa_map.txt"
   File map_vcf_updog = "~{methodName}_vcf_updog_map.txt"
   File map_bam_polyrad = "~{methodName}_bam_polyrad_map.txt"
   File map_bam_supermassa = "~{methodName}_bam_supermassa_map.txt"
   File map_bam_updog = "~{methodName}_bam_updog_map.txt"
   File map_gusmap = "~{methodName}_vcf_gusmap_map.txt"
   File map_bam_gusmap = "~{methodName}_bam_gusmap_map.txt"
 }
}

task JointDatas{
    input{
       Array[File] depths_GQ_vcf_rds
       Array[File] depths_df_vcf_rds
       Array[File] depths_updog_vcf_rds
       Array[File] depths_supermassa_vcf_rds
       Array[File] depths_polyrad_vcf_rds
       Array[File] depths_updog_bam_rds
       Array[File] depths_supermassa_bam_rds
       Array[File] depths_polyrad_bam_rds
       Array[File] times_vcf_df
       Array[File] times_vcf_GQ
       Array[File] times_vcf_updog
       Array[File] times_vcf_supermassa
       Array[File] times_vcf_polyrad
       Array[File] times_vcf_gusmap
       Array[File] times_bam_gusmap
       Array[File] times_bam_updog
       Array[File] times_bam_supermassa
       Array[File] times_bam_polyrad
       Array[File] filters_vcf_dfAndGQ
       Array[File] filters_vcf_polyrad
       Array[File] filters_vcf_supermassa
       Array[File] filters_vcf_updog
       Array[File] filters_bam_polyrad
       Array[File] filters_bam_supermassa
       Array[File] filters_bam_updog
       Array[File] RData_vcf_df
       Array[File] RData_vcf_GQ
       Array[File] RData_vcf_polyrad
       Array[File] RData_vcf_supermassa
       Array[File] RData_vcf_updog
       Array[File] RData_bam_polyrad
       Array[File] RData_bam_supermassa
       Array[File] RData_bam_updog
       Array[File] RData_gusmap
       Array[File] RData_bam_gusmap
       Array[File] map_vcf_df
       Array[File] map_vcf_GQ
       Array[File] map_vcf_polyrad
       Array[File] map_vcf_supermassa
       Array[File] map_vcf_updog
       Array[File] map_bam_polyrad
       Array[File] map_bam_supermassa
       Array[File] map_bam_updog
       Array[File] map_gusmap
       Array[File] map_bam_gusmap
    }

    command <<<

       R --vanilla --no-save <<RSCRIPT

       system("ln -s ~{sep = " . ; ln -s " depths_GQ_vcf_rds} .")
       system("ln -s ~{sep = " . ; ln -s " depths_df_vcf_rds} .")
       system("ln -s ~{sep = " . ; ln -s " depths_updog_vcf_rds } .")
       system("ln -s ~{sep = " . ; ln -s " depths_supermassa_vcf_rds } .")
       system("ln -s ~{sep = " . ; ln -s " depths_polyrad_vcf_rds } .")
       system("ln -s ~{sep = " . ; ln -s " depths_updog_bam_rds } .")
       system("ln -s ~{sep = " . ; ln -s " depths_supermassa_bam_rds } .")
       system("ln -s ~{sep = " . ; ln -s " depths_polyrad_bam_rds } .")
       system("ln -s ~{sep = " . ; ln -s " times_vcf_df } .")
       system("ln -s ~{sep = " . ; ln -s " times_vcf_GQ } .")
       system("ln -s ~{sep = " . ; ln -s " times_vcf_updog } .")
       system("ln -s ~{sep = " . ; ln -s " times_vcf_supermassa } .")
       system("ln -s ~{sep = " . ; ln -s " times_vcf_polyrad } .")
       system("ln -s ~{sep = " . ; ln -s " times_vcf_gusmap } .")
       system("ln -s ~{sep = " . ; ln -s " times_bam_updog } .")
       system("ln -s ~{sep = " . ; ln -s " times_bam_supermassa } .")
       system("ln -s ~{sep = " . ; ln -s " times_bam_polyrad } .")
       system("ln -s ~{sep = " . ; ln -s " times_bam_gusmap } .")
       system("ln -s ~{sep = " . ; ln -s " filters_vcf_dfAndGQ } .")
       system("ln -s ~{sep = " . ; ln -s " filters_vcf_polyrad } .")
       system("ln -s ~{sep = " . ; ln -s " filters_vcf_supermassa } .")
       system("ln -s ~{sep = " . ; ln -s " filters_vcf_updog } .")
       system("ln -s ~{sep = " . ; ln -s " filters_bam_polyrad } .")
       system("ln -s ~{sep = " . ; ln -s " filters_bam_supermassa } .")
       system("ln -s ~{sep = " . ; ln -s " filters_bam_updog } .")
       system("ln -s ~{sep = " . ; ln -s " RData_vcf_df } .")
       system("ln -s ~{sep = " . ; ln -s " RData_vcf_GQ } .")
       system("ln -s ~{sep = " . ; ln -s " RData_vcf_polyrad } .")
       system("ln -s ~{sep = " . ; ln -s " RData_vcf_supermassa } .")
       system("ln -s ~{sep = " . ; ln -s " RData_vcf_updog } .")
       system("ln -s ~{sep = " . ; ln -s " RData_bam_polyrad } .")
       system("ln -s ~{sep = " . ; ln -s " RData_bam_supermassa } .")
       system("ln -s ~{sep = " . ; ln -s " RData_bam_updog } .")
       system("ln -s ~{sep = " . ; ln -s " RData_gusmap } .")
       system("ln -s ~{sep = " . ; ln -s " RData_bam_gusmap } .")
       system("ln -s ~{sep = " . ; ln -s " map_vcf_df } .")
       system("ln -s ~{sep = " . ; ln -s " map_vcf_GQ } .")
       system("ln -s ~{sep = " . ; ln -s " map_vcf_polyrad } .")
       system("ln -s ~{sep = " . ; ln -s " map_vcf_supermassa } .")
       system("ln -s ~{sep = " . ; ln -s " map_vcf_updog } .")
       system("ln -s ~{sep = " . ; ln -s " map_bam_polyrad } .")
       system("ln -s ~{sep = " . ; ln -s " map_bam_supermassa } .")
       system("ln -s ~{sep = " . ; ln -s " map_bam_updog } .")
       system("ln -s ~{sep = " . ; ln -s " map_gusmap } .")
       system("ln -s ~{sep = " . ; ln -s "  map_bam_gusmap } .")

       ## Map, times and RDatas

       SNPCall <- c("gatk", "freebayes")
       GenoCall <- c("supermassa", "updog", "df", "polyrad", "GQ", "gusmap")
       CountsFrom <- c("bam", "vcf")

       data_times <- data_map <- temp <- temp.1 <- data.frame()
       z <- 1

       names_list <- vector()
       data_RDatas <- list()
         for(i in SNPCall){
            for(j in GenoCall){
               for(w in CountsFrom){
                 if((j == "df" | j =="GQ") && w=="bam"){
                 } else {
                  temp <- read.table(paste0(i,"_",w,"_",j,"_map.txt"), header=T)
                  temp.1 <- read.table(paste0(i,"_",w,"_",j,"_times.txt"), header=T)
                  load(paste0(i,"_",w,"_",j,".RData"))
                  data_map <- rbind(data_map, temp)
                  data_times <- rbind(data_times, temp.1)
                  data_RDatas[[z]] <- map_out
                  names_list <- c(names_list,paste0(i,"_",w,"_",j))
                  z <- z + 1
                 }
                }
             }
          }

        saveRDS(data_map, "data_map.rds")
        saveRDS(data_times, "data_times.rds")

        names(data_RDatas) <- names_list
        save(data_RDatas, file = "data_RDatas.RData")

        # Filters

        GenoCall <- c("supermassa", "updog", "dfAndGQ", "polyrad")

        data_filters <- temp <- data.frame()
        for(i in SNPCall){
            for(j in GenoCall){
                for(w in CountsFrom){
                    if((j == "dfAndGQ") & w=="bam"){
                    } else {
                     temp <- read.table(paste0(i,"_filters_",w,"_",j,".txt"), header=T)
                    }
                  data_filters <- rbind(data_filters, temp)
                }
            }
        }

        saveRDS(data_filters, "data_filters.rds")

       # Depths

        GenoCall <- c("supermassa", "updog", "df", "polyrad", "GQ")

        data_depths <- data.frame()
        for(i in SNPCall){
            for(j in GenoCall){
                for(w in CountsFrom){
                   if((j == "df" | j == "GQ") & w == "bam"){
                   } else {
                     temp <- readRDS(paste0(i,"_", j, "_", w, "_depths.rds"))
                     temp <- data.frame(SNPCall=i, GenoCall=j, CountsFrom=w, temp)
                     data_depths <- rbind(data_depths, temp)
                   }
                }
            }
        }

       saveRDS(data_depths, "data_depths.rds")

       RSCRIPT
     >>>

    runtime{
      docker:"taniguti/onemap"
    }

    output{
      File data_depths = "data_depths.rds"
      File data_filters = "data_filters.rds"
      File data_map = "data_map.rds"
      File data_times = "data_times.rds"
      File data_RDatas = "data_RDatas.RData"
    }
}
