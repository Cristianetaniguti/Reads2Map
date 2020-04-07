version 1.0

import "./utilsR.wdl" as utilsR

workflow GusmapMaps{
  input{
    File simu_onemap_obj
    File vcf_file
    File new_vcf_file
    String SNPCall_program
    String GenotypeCall_program
    File tot_mks
    File real_phases
    String cMbyMb
  }

  Array[String] counts                      = ["vcf", "bam"]
  Array[File] vcfs                          = [vcf_file, new_vcf_file]
  Array[Pair[String, File]] counts_and_vcfs = zip(counts, vcfs)
  
  scatter(vcf in counts_and_vcfs){
    call GusmapReport{
        input:
          vcf_file = vcf.right,
          simu_onemap_obj = simu_onemap_obj, 
          SNPCall_program = SNPCall_program,
          GenotypeCall_program = GenotypeCall_program,
          CountsFrom = vcf.left,
          tot_mks = tot_mks,
          real_phases = real_phases,
          cMbyMb = cMbyMb
        }
    }
     
   output{
      Array[File] RDatas = GusmapReport.maps_RData
      Array[File] maps_report = GusmapReport.maps_report
      Array[File] times = GusmapReport.times
   }
}

task GusmapReport{
  input{
    File vcf_file
    File simu_onemap_obj
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    File tot_mks
    File real_phases
    String cMbyMb
  }
  
  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)
      library(GUSMap)
      source("/opt/scripts/functions_simu.R")
      
      simu_onemap_obj <- load("~{simu_onemap_obj}")
      simu_onemap_obj <- get(simu_onemap_obj)
      
      if(tail(strsplit("~{vcf_file}", "[.]")[[1]],1) =="gz") {
          vcf.temp <- paste0("~{SNPCall_program}",".vcf")
          system(paste0("zcat ", "~{vcf_file}", " > ", vcf.temp))
          vcf_file <- vcf.temp
       } else {
          vcf_file <- "~{vcf_file}"
       }
      
      tot_mks <- read.table("~{tot_mks}")
      real_phases <- read.table("~{real_phases}")
      
      times <- system.time(create_gusmap_report(vcf_file, gab= simu_onemap_obj,"~{SNPCall_program}", 
                                                     "~{GenotypeCall_program}", TRUE, "~{CountsFrom}", tot_mks,real_phases, ~{cMbyMb}))

      outname <- paste0("map_", "~{SNPCall_program}", "_", "~{CountsFrom}", "_", "~{GenotypeCall_program}", "_", TRUE)
      times <- data.frame(meth = outname, time = times[3])
      
      times_temp <- system.time(create_gusmap_report(vcf_file, gab= simu_onemap_obj,"~{SNPCall_program}", 
                                                     "~{GenotypeCall_program}", FALSE, "~{CountsFrom}", tot_mks,real_phases, ~{cMbyMb}))

      outname <- paste0("map_", "~{SNPCall_program}", "_", "~{CountsFrom}", "_", "~{GenotypeCall_program}", "_", FALSE)
      times_temp <- data.frame(meth = outname, time = times_temp[3])
      
    # Joint maps data.frames
      map_temp <- read.table("map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_TRUE.txt")
      map_joint <- map_temp
      map_temp <- read.table("map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_FALSE.txt")
      map_joint <- rbind(map_joint, map_temp)
      write.table(map_joint, file = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt", col.names = F)
      
      # Joint RDatas
      RDatas_joint <- list()
      map_temp <- load("map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_TRUE.RData")
      RDatas_joint[[1]] <- get(map_temp)
      map_temp <- load("map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_FALSE.RData")
      RDatas_joint[[2]] <- get(map_temp)
      names(RDatas_joint) <- c("map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_TRUE", "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_FALSE")
      save(RDatas_joint, file= "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.RData")
      
      # Joint times data.frames
      times <- rbind(times, times_temp)
      write.table(times, "times_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt", col.names = F)
  
    RSCRIPT
  >>>
  
  runtime{
    docker: "taniguti/onemap"
    time:"72:00:00"
    mem:"--nodes=1"
    cpu:1
  }
  
  output{
    File maps_report = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt" 
    File maps_RData = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.RData"
    File times = "times_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
  }

}
