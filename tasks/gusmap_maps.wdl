version 1.0


workflow GusmapMaps {
  input {
    File simu_onemap_obj
    File vcf_file
    File new_vcf_file
    String SNPCall_program
    String GenotypeCall_program
    File ref_alt_alleles
    File simulated_phases
  }

  Array[String] counts = ["vcf", "bam"]
  Array[File] vcfs = [vcf_file, new_vcf_file]
  Array[Pair[String, File]] counts_and_vcfs = zip(counts, vcfs)

  scatter (vcf in counts_and_vcfs) {
    call GusmapReport {
        input:
          vcf_file = vcf.right,
          simu_onemap_obj = simu_onemap_obj,
          SNPCall_program = SNPCall_program,
          GenotypeCall_program = GenotypeCall_program,
          CountsFrom = vcf.left,
          ref_alt_alleles = ref_alt_alleles,
          simulated_phases = simulated_phases
        }
    }

   output {
      Array[File] RDatas = GusmapReport.maps_RData
      Array[File] maps_report = GusmapReport.maps_report
      Array[File] times = GusmapReport.times
   }
}

task GusmapReport {
  input {
    File vcf_file
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    File simu_onemap_obj
    File ref_alt_alleles
    File simulated_phases
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

      ref_alt_alleles <- read.table("~{ref_alt_alleles}")
      simulated_phases <- read.table("~{simulated_phases}")

      times_fake <- system.time(create_gusmap_report(vcf_file, gab= simu_onemap_obj,"~{SNPCall_program}",
                                                     "~{GenotypeCall_program}", TRUE, "~{CountsFrom}", ref_alt_alleles,simulated_phases))

      outname <- paste0("map_", "~{SNPCall_program}", "_", "~{CountsFrom}", "_", "~{GenotypeCall_program}", "_", TRUE)
      times <- data.frame(meth = outname, time = times[3])

      # If there is no fake, map will not run again
      if(all(info_fake[[2]][,8])){
        info_correct <- info_fake
        times_temp <- times_fake
        
        est.pos <- info_fake[[2]][,2]
        real.type <- rep(NA, nrow(info_correct[[2]]))
        temp.type <- simu_onemap_obj$segr.type[which(simu_onemap_obj$POS %in% est.pos)]
        real.type[which(est.pos %in% as.character(simu_onemap_obj$POS))] <- temp.type
        real.type[which(is.na(real.type))] <- "non-informative"
        poscM <- ref_alt_alleles$pos.map[which(as.numeric(as.character(ref_alt_alleles$pos)) %in% as.numeric(as.character(est.pos)))]
        poscM.norm <- poscM-poscM[1]
        diff <- sqrt((poscM.norm - info_fake[[2]][,3])^2)
        real.phase <- simulated_phases[which(simulated_phases$pos%in%est.pos),][,2]
        
        info_correct[[2]]$real.type <- real.type
        info_correct[[2]]$real.phases <- real.phase
        info_correct[[2]]$poscM <- poscM
        info_correct[[2]]$poscM.norm <- poscM.norm
        info_correct[[2]]$diff <- diff
        
      } else {
        times_temp <- system.time(info_correct <- create_gusmap_report(vcf_file, gab= simu_onemap_obj,"gatk",
                                                      "gusmap", FALSE, "vcf", ref_alt_alleles,simulated_phases))
        
        times_temp <- data.frame(meth = outname, time = times_temp[3])
      }

      outname <- paste0("map_", "~{SNPCall_program}", "_", "~{CountsFrom}", "_", "~{GenotypeCall_program}", "_", FALSE)

      # Joint maps data.frames
      map_joint <- rbind(info_fake[[2]], info_correct[[2]])
            write.table(map_joint, file = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt", col.names = F)

      # Joint RDatas
      RDatas_joint <- list()
      RDatas_joint[[1]] <- info_fake[[1]]
      RDatas_joint[[2]] <- info_correct[[1]]

      names(RDatas_joint) <- c("map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_TRUE", 
                              "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_FALSE")
      save(RDatas_joint, file= "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.RData")

      # Joint times data.frames
      times <- rbind(times, times_temp)
      write.table(times, "times_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt", col.names = F)

    RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/onemap_workflows"
    preemptible: 3
    memory: "4 GB"
    cpu:1
  }

  output {
    File maps_report = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
    File maps_RData = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.RData"
    File times = "times_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
  }
}
