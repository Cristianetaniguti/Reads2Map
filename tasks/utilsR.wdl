version 1.0

task vcf2onemap{
   input{
     File vcf_file
     String cross
     String SNPCall_program
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
          save(onemap.obj, file=paste0("~{SNPCall_program}", "_vcf", "_onemap.obj.RData"))
                
        RSCRIPT
    >>>
    runtime{
      docker:"taniguti/onemap"
      time:"10:00:00"
      mem:"--nodes=1"
      cpu:1
    }

    output{
      File onemap_obj = "~{SNPCall_program}_vcf_onemap.obj.RData"
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
    time:"10:00:00"
    mem:"--nodes=1"
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
    docker: "taniguti/onemap"
    time:"72:00:00"
    mem:"--nodes=1"
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
    time:"24:00:00"
    mem:"--nodes=1"
    cpu:1
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
    time:"10:00:00"
    mem:"--nodes=1"
    cpu:1
  }
  
  output{
    File error_onemap_obj = "onemap_obj_globalError.RData"
  }
}


task BamDepths2Vcf{
  input{
    String SNPCall_program
    File vcf_file
    File freebayes_ref_bam 
    File freebayes_alt_bam 
    File gatk_ref_bam 
    File gatk_alt_bam 
    File gatk_example_alleles
    File freebayes_example_alleles
  }
  
  command <<<
    R --vanilla --no-save <<RSCRIPT
    
      library(onemap)
      library(vcfR)
      source("/opt/scripts/functions_simu.R")
       
      system("cp ~{freebayes_ref_bam} .")
      system("cp ~{freebayes_alt_bam} .")
      system("cp ~{gatk_ref_bam} .")
      system("cp ~{gatk_alt_bam} .")
      system("cp ~{gatk_example_alleles} .")
      system("cp ~{freebayes_example_alleles} .")
      
       ## Depths from bam
       depths.alt <- read.table(paste0("~{SNPCall_program}", "_alt_depth_bam.txt"), header = T)
       depths.ref <- read.table(paste0("~{SNPCall_program}", "_ref_depth_bam.txt"), header = T)

       depths <- list("ref" = depths.ref, "alt"=depths.alt)

       if(tail(strsplit("~{vcf_file}", "[.]")[[1]],1) =="gz") {
          vcf.temp <- paste0("~{SNPCall_program}",".", sample(1000,1), ".vcf")
          system(paste0("zcat ", "~{vcf_file}", " > ", vcf.temp))
          vcf_file <- vcf.temp
       }
      
       allele_file <- paste0("~{SNPCall_program}","_example4ref_alt_alleles.txt")
       bam_vcf <- make_vcf(vcf_file, depths, "~{SNPCall_program}", allele_file, "bam_vcf.vcf")
       
       bam_vcfR <- read.vcfR(bam_vcf)
       save(bam_vcfR, file="bam_vcfR.RData")
       

    RSCRIPT
  >>>
  
  runtime{
    docker:"taniguti/onemap"
    time:"15:00:00"
    mem:"--nodes=1"
    cpu:1
  }
  
  output{
    File bam_vcf = "bam_vcf.vcf"
    File bam_vcfR = "bam_vcfR.RData"
  }
}








