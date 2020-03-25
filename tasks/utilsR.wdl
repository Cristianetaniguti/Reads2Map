version 1.0

task vcf2onemap{
   input{
     File vcf_file
     String cross
     String SNPCall_program
     String GenotypeCall_program
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
                                 parent1="P1",
                                 parent2="P2",
                                 f1 = f1)
                save(onemap.obj, file=paste0("~{SNPCall_program}", "_", 
                "~{GenotypeCall_program}", "_onemap.obj.RData"))
                
        RSCRIPT
    >>>
    runtime{
      docker:"taniguti/onemap"
    }

    output{
      File onemap_obj = "~{SNPCall_program}_~{GenotypeCall_program}_onemap.obj.RData"
      File vcfR_obj = "vcfR_obj.RData"
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
    docker: "taniguti/onemap"
  }
  
  output{
    File filters_report = "filters_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
    File onemap_obj_filtered = "onemap_obj_filtered.RData"
  }
}

task MapsReport{
  input{
   File onemap_obj
   File tot_mks
   File simu_onemap_obj
   String SNPCall_program
   String GenotypeCall_program
   String CountsFrom
   String cMbyMb
   File real_phases
  }
  
  command <<<
      R --vanilla --no-save <<RSCRIPT
      library(onemap)
      source("/opt/scripts/functions_simu.R")
      
      filtered_onemap <- load("~{onemap_obj}")
      filtered_onemap <- get(filtered_onemap)
      
      simu_onemap_obj <- load("~{simu_onemap_obj}")
      simu_onemap_obj <- get(simu_onemap_obj)
      tot_mks <- read.table("~{tot_mks}")
      real_phases <- read.table("~{real_phases}")

      ## Without false SNPs
      times <-system.time(create_maps_report(input.seq = filtered_onemap,
                                             tot_mks = tot_mks, gab = simu_onemap_obj,
                                             "~{SNPCall_program}" , "~{GenotypeCall_program}",
                                             fake= F, "~{CountsFrom}", ~{cMbyMb}, real_phases))

      outname <- paste0("map_", "~{SNPCall_program}", "_", "~{CountsFrom}", 
                        "_", "~{GenotypeCall_program}", "_", FALSE)
      
      times <- data.frame(meth = outname, time = times[3])

      
      ## With false SNPs
      times_temp <-system.time(create_maps_report(input.seq = filtered_onemap,
                                             tot_mks = tot_mks, gab = simu_onemap_obj,
                                             "~{SNPCall_program}" , "~{GenotypeCall_program}",
                                             fake= T, "~{CountsFrom}", ~{cMbyMb}, real_phases))

      outname <- paste0("map_", "~{SNPCall_program}", "_", "~{CountsFrom}", 
                        "_", "~{GenotypeCall_program}", "_", TRUE)
                        
      # Joint maps data.frames
      map_temp <- read.table("map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_TRUE.txt")
      map_joint <- map_temp
      map_temp <- read.table("map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_FALSE.txt")
      map_joint <- rbind(map_joint, map_temp)
      write.table(map_joint, file = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt")
      
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
      write.table(times, "times.txt")
      
      RSCRIPT
  >>>
  
  runtime{
    docker: "taniguti/onemap" 
  }
  
  output{
    File maps_report = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt" 
    File maps_RData = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.RData"
    File times = "times.txt"
  }
}

task ErrorsReport{
  input{
    File onemap_obj
    File simu_onemap_obj
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
  }
  
  command <<<
      R --vanilla --no-save <<RSCRIPT
        library(onemap)
        source("/opt/scripts/functions_simu.R")
        
        onemap_obj <- load("~{onemap_obj}")
        onemap_obj <- get(onemap_obj)
        
        simu_onemap_obj <- load("~{simu_onemap_obj}")
        simu_onemap_obj <- get(simu_onemap_obj)
        
        create_errors_report(onemap_obj = onemap_obj, simu_onemap_obj,
                             "~{SNPCall_program}" , "~{GenotypeCall_program}",
                             "~{CountsFrom}")
  
      RSCRIPT
  >>>
  
  runtime{
    docker: "taniguti/onemap"
  }
  
  output{
    File errors_report = "errors_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
  }
}

task GlobalError{
  input{
    File onemap_obj 
    String SNPCall_program 
    String GenotypeCall_program 
    String CountsFrom 
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
  runtime{
    docker: "taniguti/onemap"
  }
  
  output{
    File error_onemap_obj = "onemap_obj_globalError.RData"
  }
}

task GQProbs{
  input{
    File vcfR_obj
    File onemap_obj
    String cross
  }
  
  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)
      library(vcfR)
      
      cross <- "~{cross}"
        
      if(cross == "F1"){
         f1 = NULL
      } else if (cross == "F2"){
         f1 = "F1"
      }
      
      vcf_temp <- load("~{vcfR_obj}")
      vcf <- get(vcf_temp)
      
      onemap_obj <- load("~{onemap_obj}")
      onemap_obj <- get(onemap_obj)
      
      # MAPS REPORT - GQ
      gq <- extract_depth(vcfR.object=vcf,
                               onemap.object=onemap_obj,
                               vcf.par="GQ",
                               parent1="P1",
                               parent2="P2",
                               f1 = f1,
                               recovering=FALSE)
          
      gq_onemap_obj <- create_probs(onemap.obj = onemap_obj, genotypes_errors=gq)
      save(gq_onemap_obj, file="gq_onemap_obj.RData")
      
    RSCRIPT
  >>>
  runtime{
    docker:"taniguti/onemap"
  }
  
  output{
    File gq_onemap_obj = "gq_onemap_obj.RData"
  }
}

task UpdogProbs{
  input{
    File vcfR_obj
    File onemap_obj
  }
  
  command <<<
     R --vanilla --no-save <<RSCRIPT
       library(onemap)
        
       vcf_temp <- load("~{vcfR_obj}")
       vcf <- get(vcf_temp)
       
       onemap_obj_temp <- load("~{onemap_obj}")
       onemap_obj <- get(onemap_obj_temp)
       
       updog_onemap_obj <- updog_genotype(vcfR.object=vcf,
                                    onemap.object=onemap_obj,
                                    vcf.par="AD",
                                    parent1="P1",
                                    parent2="P2",
                                    f1 = f1,
                                    recovering=TRUE,
                                    mean_phred=20,
                                    cores=3,
                                    depths=NULL,
                                    global_error = NULL,
                                    use_genotypes_errors = FALSE,
                                    use_genotypes_probs = TRUE)
       
       save(updog_onemap_obj, file="updog_onemap_obj.RData")
  
     RSCRIPT
  >>>
  
  runtime{
    docker:"taniguti/onemap"
  }
  
  output{
    File updog_onemap_obj = "updog_onemap_obj.RData"
  }

}

task SupermassaProbs{
  input{
    File vcfR_obj
    File onemap_obj
  }
  
  command <<<
     R --vanilla --no-save <<RSCRIPT
       library(onemap)
        
       vcf_temp <- load("~{vcfR_obj}")
       vcf <- get(vcf_temp)
       
       onemap_obj_temp <- load("~{onemap_obj}")
       onemap_obj <- get(onemap_obj_temp)
       
       supermassa_onemap_obj <- supermassa_genotype(vcfR.object=vcf,
                                    onemap.object=onemap_obj,
                                    vcf.par="AD",
                                    parent1="P1",
                                    parent2="P2",
                                    f1 = f1,
                                    recovering=TRUE,
                                    mean_phred=20,
                                    cores=3,
                                    depths=NULL,
                                    global_error = NULL,
                                    use_genotypes_errors = FALSE,
                                    use_genotypes_probs = TRUE)
       
       save(supermassa_onemap_obj, file="supermassa_onemap_obj.RData")
  
     RSCRIPT
  >>>
  
  runtime{
    docker:"taniguti/onemap"
  }
  
  output{
    File supermassa_onemap_obj = "supermassa_onemap_obj.RData"
  }

}

task PolyradProbs{
  input{
    File vcf_file
    File onemap_obj
    String cross
  }
  
  command <<<
     R --vanilla --no-save <<RSCRIPT
       library(onemap)
        
       onemap_obj_temp <- load("~{onemap_obj}")
       onemap_obj <- get(onemap_obj_temp)

       polyrad_onemap_obj <- polyRAD_genotype(vcf="~{vcf_file}",
                                        onemap.obj=onemap_obj,
                                        parent1="P1",
                                        parent2="P2",
                                        f1 = f1,
                                        crosstype="~{cross}",
                                        global_error = NULL,
                                        use_genotypes_errors = FALSE,
                                        use_genotypes_probs = TRUE)
       
       save(polyrad_onemap_obj, file="polyrad_onemap_obj.RData")
  
     RSCRIPT
  >>>
  
  runtime{
    docker:"taniguti/onemap"
  }
  
  output{
    File polyrad_onemap_obj = "polyrad_onemap_obj.RData"
  }

}


