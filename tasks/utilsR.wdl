version 1.0

task vcf2onemap{
   input{
     File vcf_file
     String cross
     String SNPCall_program
     String parent1
     String parent2
   }

   command <<<

        R --vanilla --no-save <<RSCRIPT
          library(onemap)
          library(vcfR)

          cross <- "~{cross}"

          if(cross == "F1"){
            cross <- "outcross"
            f1 = NULL
          } else if (cross == "F2"){
            cross <- "f2 intercross"
            f1 = "F1"
          }

          ## READING VCF FROM PIPELINE
          vcf <- read.vcfR("~{vcf_file}")
          save(vcf, file="vcfR_obj.RData")

          onemap.obj <- onemap_read_vcfR(vcfR.object=vcf,
                                 cross= cross,
                                 parent1="~{parent1}",
                                 parent2="~{parent2}",
                                 f1 = f1)
                                 
          save(onemap.obj, file=paste0("~{SNPCall_program}", "_vcf", "_onemap.obj.RData"))

        RSCRIPT

    >>>
    runtime{
      docker:"cristaniguti/onemap_workflows"
      time:"10:00:00"
      mem:"30GB"
      cpu:1
    }

    output{
      File onemap_obj = "~{SNPCall_program}_vcf_onemap.obj.RData"
      File vcfR_obj = "vcfR_obj.RData"
    }
}


task MultiVcf2onemap{
   input{
     File? multi
     String cross
     String SNPCall_program
     String parent1
     String parent2
     Int? seed
     Int? depth 
   }

   command <<<

        R --vanilla --no-save <<RSCRIPT
          library(onemap)
          library(vcfR)

          cross <- "~{cross}"

          if(cross == "F1"){
            cross <- "outcross"
            f1 = NULL
          } else if (cross == "F2"){
            cross <- "f2 intercross"
            f1 = "F1"
          }

          vcf_file <- "~{multi}"
          ## READING VCF FROM PIPELINE
          vcf <- read.vcfR(vcf_file)

          onemap.obj <- onemap_read_vcfR(vcfR.object=vcf,
                                 cross= cross,
                                 parent1="~{parent1}",
                                 parent2="~{parent2}",
                                 f1 = f1,
                                 only_biallelic = F)

          multi_names <- list()
          if(dim(onemap.obj[[1]])[2] == 0){
            multi_names[[1]] <- 0
          } else { 
            multi_names[[1]] <- colnames(onemap.obj[[1]])
          }
          names(multi_names) <- "~{depth}_~{seed}_~{SNPCall_program}"
          save(multi_names, file = "~{depth}_~{seed}_~{SNPCall_program}_multi.names.RData")  

          save(onemap.obj, file=paste0("~{SNPCall_program}", "_vcf_multi_onemap.obj.RData"))

        RSCRIPT

    >>>
    runtime{
      docker:"cristaniguti/onemap_workflows"
      time:"15:00:00"
      mem:"30GB"
      cpu:1
    }

    output{
      File onemap_obj = "~{SNPCall_program}_vcf_multi_onemap.obj.RData"
      File multi_names = "~{depth}_~{seed}_~{SNPCall_program}_multi.names.RData"
    }
}


task FiltersReport{
  input{
    File onemap_obj
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)
      source("/opt/scripts/functions_simu.R")

      temp <- load("~{onemap_obj}")
      temp.obj <- get(temp)
      onemap_obj_filtered <- create_filters_report(temp.obj, "~{SNPCall_program}",
                                           "~{CountsFrom}", "~{GenotypeCall_program}")
      save(onemap_obj_filtered, file="onemap_obj_filtered.RData")

    RSCRIPT

  >>>

  runtime{
    docker: "cristaniguti/onemap_workflows"
    time:"10:00:00"
    mem:"30GB"
    cpu:1
  }

  output{
    File filters_report = "filters_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
    File onemap_obj_filtered = "onemap_obj_filtered.RData"
  }
}

task FiltersReportEmp{
  input{
    File onemap_obj
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    String chromosome
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)
      source("/opt/scripts/functions_empirical.R")

      temp <- load("~{onemap_obj}")
      temp.obj <- get(temp)
      onemap_obj_filtered <- create_filters_report(temp.obj, "~{SNPCall_program}",
                                           "~{CountsFrom}", "~{GenotypeCall_program}", "~{chromosome}")
      save(onemap_obj_filtered, file="onemap_obj_filtered.RData")

    RSCRIPT

  >>>

  runtime{
    docker: "cristaniguti/onemap_workflows"
    time:"10:00:00"
    mem:"30GB"
    cpu:1
  }

  output{
    File filters_report = "filters_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
    File onemap_obj_filtered = "onemap_obj_filtered.RData"
  }
}




task MapsReport{
  input{
   File onemap_obj
   File ref_alt_alleles
   File simu_onemap_obj
   String SNPCall_program
   String GenotypeCall_program
   String CountsFrom
   File simulated_phases
  }

  command <<<
      R --vanilla --no-save <<RSCRIPT
      library(onemap)
      source("/opt/scripts/functions_simu.R")

      filtered_onemap <- load("~{onemap_obj}")
      filtered_onemap <- get(filtered_onemap)

      simu_onemap_obj <- load("~{simu_onemap_obj}")
      simu_onemap_obj <- get(simu_onemap_obj)
      ref_alt_alleles <- read.table("~{ref_alt_alleles}")
      simulated_phases <- read.table("~{simulated_phases}")

      ## Without false SNPs
      times <-system.time(create_maps_report(input.seq = filtered_onemap,
                                             tot_mks = ref_alt_alleles, gab = simu_onemap_obj,
                                             "~{SNPCall_program}" , "~{GenotypeCall_program}",
                                             fake= F, "~{CountsFrom}", simulated_phases))

      outname <- paste0("map_", "~{SNPCall_program}", "_", "~{CountsFrom}",
                        "_", "~{GenotypeCall_program}", "_", FALSE)

      times <- data.frame(meth = outname, time = times[3])


      ## With false SNPs
      times_temp <-system.time(create_maps_report(input.seq = filtered_onemap,
                                             tot_mks = ref_alt_alleles, gab = simu_onemap_obj,
                                             "~{SNPCall_program}" , "~{GenotypeCall_program}",
                                             fake= T, "~{CountsFrom}", simulated_phases))

      outname <- paste0("map_", "~{SNPCall_program}", "_", "~{CountsFrom}",
                        "_", "~{GenotypeCall_program}", "_", TRUE)

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
      times_temp <- data.frame(meth = outname, time = times_temp[3])
      times <- rbind(times, times_temp)
      write.table(times, "times_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt", col.names = F)

      RSCRIPT

  >>>

  runtime{
    docker: "cristaniguti/onemap_workflows"
    time:"24:00:00"
    mem:"50GB"
    cpu:4
  }

  output{
    File maps_report = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
    File maps_RData = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.RData"
    File times = "times_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
  }
}

task ErrorsReport{
  input{
    File onemap_obj
    File vcfR_obj
    File simu_onemap_obj
    File simu_vcfR
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
  }

  command <<<
      R --vanilla --no-save <<RSCRIPT
        
        library(onemap)

        temp <- load("~{onemap_obj}")
        df <- get(temp)

        temp <- load("~{vcfR_obj}")
        vcf <- get(temp)

        p <- create_depths_profile(onemap.obj = df, vcfR.object = vcf, parent1 = "P1",
        parent2 = "P2", vcf.par = "AD",recovering = FALSE, GTfrom = "vcf", alpha=0.1,
        rds.file = paste0("~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_vcf_depths.rds"))

        df <- readRDS(paste0("~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_vcf_depths.rds"))
        df <- cbind(SNPCall = "~{SNPCall_program}", CountsFrom = "~{CountsFrom}",
                    GenoCall="~{GenotypeCall_program}", df)

        simu <- load("~{simu_vcfR}")

        write.table(simu, file="errors_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt", row.names=F, quote=F, col.names=F)

        # library(onemap)
        # source("/opt/scripts/functions_simu.R")

        # onemap_obj <- load("~{onemap_obj}")
        # onemap_obj <- get(onemap_obj)

        # simu_onemap_obj <- load("~{simu_onemap_obj}")
        # simu_onemap_obj <- get(simu_onemap_obj)

        # create_errors_report(onemap_obj = onemap_obj, simu_onemap_obj,
        #                      "~{SNPCall_program}" , "~{GenotypeCall_program}",
        #                      "~{CountsFrom}")

      RSCRIPT

  >>>

  runtime{
    docker: "cristaniguti/onemap_workflows"
    time:"10:00:00"
    mem:"30GB"
    cpu:1
  }

  output{
    File errors_report = "errors_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
  }
}

task GlobalError {
  input {
    File onemap_obj
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
    library(onemap)

    onemap_obj <- load("~{onemap_obj}")
    onemap_obj <- get(onemap_obj)

    onemap_obj_globalError <- create_probs(onemap.obj = onemap_obj, global_error = 0.05)
    save(onemap_obj_globalError, file = "onemap_obj_globalError.RData")

    RSCRIPT

  >>>
  runtime {
    docker: "cristaniguti/onemap_workflows"
    time:"10:00:00"
    mem:"30GB"
    cpu:1
  }

  output{
    File error_onemap_obj = "onemap_obj_globalError.RData"
  }
}


task BamDepths2Vcf{
  input{
    File vcf_file
    File ref_bam
    File alt_bam
    File example_alleles
    String program
    Int max_cores
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT

      library(onemap)
      library(vcfR)
      library(doParallel)
      source("/opt/scripts/functions_simu.R")

      system("cp ~{ref_bam} .")
      system("cp ~{alt_bam} .") 
      system("cp ~{example_alleles} .")

       ## Depths from bam
       depths.alt <- read.table("~{alt_bam}", header = T)
       depths.ref <- read.table("~{ref_bam}", header = T)

       depths <- list("ref" = depths.ref, "alt"=depths.alt)

       if(tail(strsplit("~{vcf_file}", "[.]")[[1]],1) =="gz") {
          vcf.temp <- paste0("vcf.temp",".", sample(1000,1), ".vcf")
          system(paste0("zcat ", "~{vcf_file}", " > ", vcf.temp))
          vcf_file <- vcf.temp
       } else {
          vcf_file <- "~{vcf_file}"
       }

       allele_file <- paste0("~{example_alleles}")
       bam_vcf <- make_vcf(vcf_file, depths, allele_file, "~{program}_bam_vcf.vcf", cores = ~{max_cores})

       bam_vcfR <- read.vcfR(bam_vcf)
       save(bam_vcfR, file="~{program}_bam_vcfR.RData")

    RSCRIPT

  >>>

  runtime{
    docker:"cristaniguti/onemap_workflows"
    time:"15:00:00"
    mem:"30GB"
    cpu:1
  }

  output{
    File bam_vcf = "~{program}_bam_vcf.vcf"
    File bam_vcfR = "~{program}_bam_vcfR.RData"
  }
}


task CheckDepths{
  input{
    File onemap_obj
    File vcfR_obj
    String parent1
    String parent2
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)

      temp <- load("~{onemap_obj}")
      df <- get(temp)

      temp <- load("~{vcfR_obj}")
      vcf <- get(temp)

      p <- create_depths_profile(onemap.obj = df, vcfR.object = vcf, parent1 = "~{parent1}",
      parent2 = "~{parent2}", vcf.par = "AD",recovering = FALSE, GTfrom = "vcf", alpha=0.1,
      rds.file = paste0("~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_vcf_depths.rds"))

      df <- readRDS(paste0("~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_vcf_depths.rds"))
      df <- cbind(SNPCall = "~{SNPCall_program}", CountsFrom = "~{CountsFrom}",
                  GenoCall="~{GenotypeCall_program}", df)
      write.table(df, file="~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_vcf_depths.txt", row.names=F, quote=F, col.names=F)

      #ggsave(filename = paste0("~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_vcf_depths.png"), p)

    RSCRIPT

  >>>

  runtime{
    docker:"cristaniguti/onemap_workflows"
    time:"10:00:00"
    mem:"60GB"
    cpu:1
  }

  output{
    File errors_report = "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_vcf_depths.txt"
  }

}

task MapsReportEmp{
  input{
   File sequence_obj
   String SNPCall_program
   String GenotypeCall_program
   String CountsFrom
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)
      source("/opt/scripts/functions_empirical.R")

      temp <- load("~{sequence_obj}")
      sequence <- get(temp)

      create_map_report(input.seq = sequence, CountsFrom = "~{CountsFrom}",
                        SNPCall = "~{SNPCall_program}", GenoCall="~{GenotypeCall_program}")

    RSCRIPT

  >>>

  runtime{
    docker:"cristaniguti/onemap_workflows"
    time:"24:00:00"
    mem:"60GB"
    cpu:4
  }

  output{
    File maps_report = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
    File maps_RData = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.RData"
    File times = "times_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
  }
}


task AddMultiallelics{
  input{
   File? onemap_obj_multi
   File onemap_obj_bi
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)
      
      temp <- load("~{onemap_obj_multi}")
      onemap_obj_multi <- get(temp)
      
      temp <- load("~{onemap_obj_bi}")
      onemap_obj_bi <- get(temp)

      onemap_obj <- create_probs(onemap_obj_multi, global_error = 0.05) # All multiallelics receives global_error = 0.05
      
      onemap_both <- combine_onemap(onemap_obj_bi, onemap_obj_multi)
      
      onemap_both <- sort_by_pos(onemap_both)
      
      save(onemap_both, file=paste0("onemap_obj_both.RData"))
                
    RSCRIPT

  >>>

  runtime{
    docker:"cristaniguti/onemap_workflows"
    time:"05:00:00"
    mem:"30GB"
    cpu:4
  }

  output{
      File onemap_obj_both = "onemap_obj_both.RData"
  }
}




