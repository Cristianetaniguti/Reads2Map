version 1.0

task vcf2onemap{
   input {
     File vcf_file
     String cross
     String SNPCall_program
     String parent1
     String parent2
     String multiallelics
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
          onemap.obj <- onemap_read_vcfR(vcf="~{vcf_file}",
                                         cross= cross,
                                         parent1="~{parent1}",
                                         parent2="~{parent2}",
                                         f1 = f1)

          save(onemap.obj, file=paste0("vcf_onemap.obj.RData"))

        RSCRIPT

    >>>
    runtime {
      docker:"cristaniguti/reads2map:0.0.1"
      # preemptible: 3
      # memory: "2 GB"
      # cpu:1
      job_name: "vcf2onemap"
      node:"--nodes=1"
      mem:"--mem=10G"
      cpu:"--ntasks=1"
      time:"10:00:00"
    }

    output{
      File onemap_obj = "vcf_onemap.obj.RData"
    }
}

task FiltersReport{
  input{
    File onemap_obj
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    Int seed
    Int depth
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)
      source("/opt/scripts/functions_simu.R")

      temp <- load("~{onemap_obj}")
      temp.obj <- get(temp)
      onemap_obj_filtered <- create_filters_report(temp.obj, "~{SNPCall_program}",
                                                  "~{CountsFrom}", "~{GenotypeCall_program}", 
                                                   ~{seed}, ~{depth})

      save(onemap_obj_filtered, file="onemap_obj_filtered.RData")

    RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/reads2map:0.0.1"
    # preemptible: 3
    # memory:"3 GB"
    # cpu:1
    job_name: "Filters"
    node:"--nodes=1"
    mem:"--mem=7G"
    cpu:"--ntasks=1"
    time:"05:00:00"
  }

  output {
    File filters_report = "filters_report.tsv.gz"
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

  runtime {
    docker: "cristaniguti/reads2map:0.0.1"
    # preemptible: 3
    # memory:"4 GB"
    # cpu:1
    job_name: "Filters"
    node:"--nodes=1"
    mem:"--mem=7G"
    cpu:"--ntasks=1"
    time:"05:00:00"
  }

  output {
    File filters_report = "~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}_filters_report.tsv.gz"
    File onemap_obj_filtered = "onemap_obj_filtered.RData"
  }
}

task MapsReport{
  input {
   File onemap_obj
   File ref_alt_alleles
   File simu_onemap_obj
   String SNPCall_program
   String GenotypeCall_program
   String CountsFrom
   File simulated_phases
   Int seed
   Int depth
   Int max_cores
  }

  command <<<
      R --vanilla --no-save <<RSCRIPT
      library(onemap)
      source("/opt/scripts/functions_simu.R")

      if(~{max_cores} > 4) cores = 4 else cores = ~{max_cores}

      filtered_onemap <- load("~{onemap_obj}")
      filtered_onemap <- get(filtered_onemap)

      simu_onemap_obj <- load("~{simu_onemap_obj}")
      simu_onemap_obj <- get(simu_onemap_obj)
      ref_alt_alleles <- read.table("~{ref_alt_alleles}")
      simulated_phases <- read.table("~{simulated_phases}")

      ## With false SNPs
      times_fake <-system.time(info_fake <- create_maps_report(input.seq = filtered_onemap,
                                                  tot_mks = ref_alt_alleles, gab = simu_onemap_obj,
                                                  "~{SNPCall_program}" , "~{GenotypeCall_program}",
                                                  fake= T, "~{CountsFrom}", simulated_phases, 
                                                  ~{seed}, ~{depth}, cores))

      times <- data.frame(seed = ~{seed}, depth = ~{depth}, SNPCall = "~{SNPCall_program}", 
                          CountsFrom = "~{CountsFrom}", GenoCall =  "~{GenotypeCall_program}", fake = "with-false",
                          time = times_fake[3])

      # It will not run if all markers are true markers
      if(all(info_fake[[2]][,"real.mks"] == "true marker")){     
        cat("skip :) \n")
        times_temp <- times_fake
        info_correct <- update_fake_info(info_fake, simu_onemap_obj, ref_alt_alleles, simulated_phases)
        
      } else {
        ## Without false SNPs
        times_temp <-system.time(info_correct <- create_maps_report(input.seq = filtered_onemap,
                                              tot_mks = ref_alt_alleles, gab = simu_onemap_obj,
                                              "~{SNPCall_program}" , "~{GenotypeCall_program}",
                                              fake= F, "~{CountsFrom}", simulated_phases, 
                                              ~{seed}, ~{depth}, cores))
      }

      # Joint maps data.frames
      map_joint <- rbind(info_fake[[2]], info_correct[[2]])
      vroom::vroom_write(map_joint, "map_report.tsv.gz", num_threads = ~{max_cores})

      # Joint RDatas
      RDatas_joint <- list()
      RDatas_joint[[1]] <- info_fake[[1]]
      RDatas_joint[[2]] <- info_correct[[1]]
      names(RDatas_joint) <- c("map_~{seed}_~{depth}_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_TRUE", 
                               "map_~{seed}_~{depth}_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_FALSE")
      save(RDatas_joint, file= "map.RData")

      # Joint times data.frames
      times_temp <- data.frame(seed = ~{seed}, depth = ~{depth}, SNPCall = "~{SNPCall_program}", 
                               CountsFrom = "~{CountsFrom}", GenoCall =  "~{GenotypeCall_program}", fake = "without-false",
                               time = times_temp[3])

      times <- rbind(times, times_temp)
      vroom::vroom_write(times, "times_report.tsv.gz", num_threads = ~{max_cores})

      RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/reads2map:0.0.1"
    # preemptible: 3
    # memory: "4 GB"
    # cpu: 4
    job_name: "MapsReport"
    node:"--nodes=1"
    mem:"--mem=20G"
    cpu:"--ntasks=4"
    time:"48:00:00"
  }

  output {
    File maps_report = "map_report.tsv.gz"
    File maps_RData = "map.RData"
    File times = "times_report.tsv.gz"
  }
}

task ErrorsReport {
  input {
    File onemap_obj
    File vcfR_obj
    File simu_onemap_obj
    File simu_vcfR
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    Int seed
    Int depth
    Int max_cores
  }

  command <<<
      R --vanilla --no-save <<RSCRIPT

        library(onemap)
        library(tidyverse)

        temp <- load("~{onemap_obj}")
        df <- get(temp)

        temp <- load("~{vcfR_obj}")
        vcf <- get(temp)

        p <- create_depths_profile(onemap.obj = df, vcfR.object = vcf, parent1 = "P1",
                                   parent2 = "P2", vcf.par = "AD",recovering = FALSE, 
                                   GTfrom = "onemap", alpha=0.1,
                                   rds.file = paste0("vcf_depths.rds"))

        df <- readRDS(paste0("vcf_depths.rds"))
        df <- cbind(seed = "~{seed}"  , depth = "~{depth}", SNPCall = "~{SNPCall_program}", CountsFrom = "~{CountsFrom}",
                    GenoCall="~{GenotypeCall_program}", df)

        simu <- load("~{simu_vcfR}")
        vcf_simu <- get(simu)

        gt.simu <- vcf_simu@gt[,-1]
        gt.simu <- as.data.frame(cbind(mks = vcf_simu@fix[,3], gt.simu))

        gt.simu <- data.frame(lapply(gt.simu, function(x) {
            x %>% gsub("[|]", "/", .) %>%
            gsub("[.]","missing", .) %>%
            gsub("0/0","homozygous-ref", .) %>%
            gsub("1/1","homozygous-alt", .) %>%
            gsub("0/1","heterozygous", .) %>%
            gsub("1/0","heterozygous", .)
        }))

        dptot <- gt.simu %>%
          pivot_longer(!mks, names_to = "ind", values_to = "gabGT") %>% inner_join(df)

        idx <- match(c("A", "AB", "BA", "B"), colnames(dptot))
        dptot <- cbind(dptot, errors = apply(dptot[,idx], 1, function(x) 1 - max(x)))

        vroom::vroom_write(dptot, "errors_report.tsv.gz", num_threads = 4)

      RSCRIPT

  >>>

  runtime{
    docker: "cristaniguti/reads2map:0.0.1"
    # preemptible: 3
    # memory: "2 GB"
    # cpu:4
    job_name: "ErrorsReport"
    node:"--nodes=1"
    mem:"--mem=5G"
    cpu:"--ntasks=1"
    time:"05:00:00"
  }

  output{
    File errors_report = "errors_report.tsv.gz"
  }
}

task GlobalError {
  input {
    File vcf_file
    File onemap_obj
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
    library(onemap)

    onemap_obj <- load("~{onemap_obj}")
    onemap_obj <- get(onemap_obj)

    onemap_obj_globalError <- create_probs(input.obj = onemap_obj, global_error = 0.05)
    save(onemap_obj_globalError, file = "onemap_obj_globalError.RData")

    RSCRIPT

  >>>
  runtime {
    docker: "cristaniguti/reads2map:0.0.1"
    # preemptible: 3
    # memory: "2 GB"
    # cpu: 1
    job_name: "GlobalError"
    node:"--nodes=1"
    mem:"--mem=20G"
    cpu:"--ntasks=1"
    time:"10:00:00"
  }

  output {
    File error_onemap_obj = "onemap_obj_globalError.RData"
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
    Int? seed
    Int? depth
    Int max_cores
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
      rds.file = paste0("vcf_depths.rds"))

      df <- readRDS(paste0("vcf_depths.rds"))
      df <- cbind(SNPCall = "~{SNPCall_program}", CountsFrom = "~{CountsFrom}",
                  GenoCall="~{GenotypeCall_program}", df ~{", seed=" + seed} ~{", depth= " + depth})

      vroom::vroom_write(df, "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_errors_report.tsv.gz", num_threads = ~{max_cores})

    RSCRIPT

  >>>

  runtime{
    docker:"cristaniguti/reads2map:0.0.1"
    # preemptible: 3
    # memory:"3 GB"
    # cpu:4
    job_name: "CheckDepths"
    node:"--nodes=1"
    mem:"--mem=20G"
    cpu:"--ntasks=1"
    time:"10:00:00"
  }

  output{
    File errors_report = "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_errors_report.tsv.gz"
  }

}

task MapsReportEmp{
  input{
   File sequence_obj
   String SNPCall_program
   String GenotypeCall_program
   String CountsFrom
   Int max_cores
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)
      source("/opt/scripts/functions_empirical.R")

      temp <- load("~{sequence_obj}")
      sequence <- get(temp)

      if(~{max_cores} > 4) cores = 4 else cores = ~{max_cores}

      times_temp <- system.time(df <- create_map_report(input.seq = sequence, CountsFrom = "~{CountsFrom}",
                                      SNPCall = "~{SNPCall_program}", GenoCall="~{GenotypeCall_program}", max_cores = cores))

      vroom::vroom_write(df[[2]], "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_map_report.tsv.gz", num_threads = ~{max_cores})
      map_out <- df[[1]]
      save(map_out,  file = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.RData")
    
      times <- data.frame(SNPCall = "~{SNPCall_program}", 
                          CountsFrom = "~{CountsFrom}", 
                          GenoCall =  "~{GenotypeCall_program}",
                          time = times_temp[3])   

      vroom::vroom_write(times, "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_times_report.tsv.gz", num_threads = ~{max_cores})

    RSCRIPT
  >>>

  runtime{
    docker:"cristaniguti/reads2map:0.0.1"
    # preemptible: 3
    # memory:"8 GB"
    # cpu:4
    job_name: "MapsEmp"
    node:"--nodes=1"
    mem:"--mem=30G"
    cpu:"--ntasks=4"
    time:"48:00:00"
  }

  output{
    File maps_report = "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_map_report.tsv.gz"
    File maps_RData = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.RData"
    File times = "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_times_report.tsv.gz"
  }
}

# Simulations taks
task JointReports{
  input {
    Array[File] default_RDatas
    Array[File] default_maps_report
    Array[File] default_filters_report
    Array[File] default_errors_report
    Array[File] default_times_report
    Array[File] SNPCaller_RDatas
    Array[File] SNPCaller_maps_report
    Array[File] SNPCaller_filters_report
    Array[File] SNPCaller_errors_report
    Array[File] SNPCaller_times_report
    Array[File] Updog_RDatas
    Array[File] Updog_maps_report
    Array[File] Updog_filters_report
    Array[File] Updog_errors_report
    Array[File] Updog_times_report
    Array[File] Polyrad_RDatas
    Array[File] Polyrad_maps_report
    Array[File] Polyrad_filters_report
    Array[File] Polyrad_errors_report
    Array[File] Polyrad_times_report
    Array[File] Supermassa_RDatas
    Array[File] Supermassa_maps_report
    Array[File] Supermassa_filters_report
    Array[File] Supermassa_errors_report
    Array[File] Supermassa_times_report
    Array[File] Gusmap_RDatas
    Array[File] Gusmap_maps_report
    Array[File] Gusmap_times_report
    Array[File] multi_names
    File Freebayes_eval
    File GATK_eval
    Int max_cores
    Int seed
    Int depth
  }

  command <<<
     R --vanilla --no-save <<RSCRIPT

      # packages
      library(tidyr)
      library(stringr)
      library(vroom)
      library(largeList)

      # function
      joint_reports <- function(default, 
                                snpcaller, 
                                updog, 
                                polyrad, 
                                supermassa,
                                gusmap = NULL){

        default    <- str_split(default, ";", simplify = T)
        SNPCaller  <- str_split(snpcaller, ";", simplify = T)
        Updog      <- str_split(updog, ";", simplify = T)
        Polyrad    <- str_split(polyrad, ";", simplify = T)
        Supermassa <- str_split(supermassa, ";", simplify = T)

        if(is.null(gusmap)){
          files <- c(default, SNPCaller, Updog, Polyrad, Supermassa)
        } else {
          Gusmap <- str_split(gusmap, ";", simplify = T)
          files <- c(default, SNPCaller, Updog, Polyrad, Supermassa, Gusmap)
        }

        joint <- vroom(files, num_threads = ~{max_cores})
        return(joint)
      }

      #########################################################################################
      # Table1: GenoCall; mks; ind; SNPCall; CountsFrom; alt; ref; gt.onemap; gt.onemap.ref.alt; 
      # gt.vcf; gt.vcf.ref.alt; gabGT; A; AB; BA; B; errors; seed; depth
      #########################################################################################

      errors_report <- joint_reports(default = "~{sep=";" default_errors_report}", 
                                     snpcaller = "~{sep=";" SNPCaller_errors_report}", 
                                     updog = "~{sep=";" Updog_errors_report}", 
                                     polyrad = "~{sep=";" Polyrad_errors_report}", 
                                     supermassa = "~{sep=";" Supermassa_errors_report}")

      ##################################################################################
      # Table2: seed; depth; CountsFrom; GenoCall; SNPCall; MK; rf; phases; real.phases;
      # real.type; real.mks; fake; poscM; poscM.norm; diff
      ##################################################################################

      maps_report <- joint_reports(default = "~{sep=";" default_maps_report}", 
                                   snpcaller = "~{sep=";" SNPCaller_maps_report}", 
                                   updog = "~{sep=";" Updog_maps_report}", 
                                   polyrad = "~{sep=";" Polyrad_maps_report}", 
                                   supermassa = "~{sep=";" Supermassa_maps_report}",
                                   gusmap = "~{sep=";" Gusmap_maps_report}")

      # Add multiallelics tag
      multi_temp <- "~{sep=";" multi_names}"
      multi_temp <- str_split(multi_temp, ";", simplify = T)

      multi_names_seed <- list()
      for(i in 1:length(multi_temp)){
        multi_temp2 <- load(multi_temp[i])
        multi_temp3 <- get(multi_temp2)
        multi_names_seed <- c(multi_names_seed, multi_temp3)
      }

      snpcall_names <- str_split(names(multi_names_seed), pattern = "_", simplify = T)
      maps_report <- as.data.frame(maps_report)

      for(i in 1:length(multi_names_seed)){
        maps_report[,"real.mks"][which(maps_report[,"seed"] == snpcall_names[i,2] & 
                                    maps_report[,"depth"] == snpcall_names[i,1] &
                                    maps_report[,"SNPCall"] == snpcall_names[i,3] &  
                                    maps_report[,"mk.name"] %in% multi_names_seed[[i]])] <- "multiallelic"
      }

      ################################################################################
      # Table3: CountsFrom; seed; depth; SNPCall; GenoCall; n_mks; distorted; 
      # redundant; mis; after
      ################################################################################

      filters_report <- joint_reports(default = "~{sep=";" default_filters_report}", 
                                   snpcaller = "~{sep=";" SNPCaller_filters_report}", 
                                   updog = "~{sep=";" Updog_filters_report}", 
                                   polyrad = "~{sep=";" Polyrad_filters_report}", 
                                   supermassa = "~{sep=";" Supermassa_filters_report}")

      #################################################################################
      # Table4: seed; depth; CountsFrom; SNPCall; GenoCall; fake; times
      #################################################################################

      times_report <- joint_reports(default = "~{sep=";" default_times_report}", 
                                   snpcaller = "~{sep=";" SNPCaller_times_report}", 
                                   updog = "~{sep=";" Updog_times_report}", 
                                   polyrad = "~{sep=";" Polyrad_times_report}", 
                                   supermassa = "~{sep=";" Supermassa_times_report}",
                                   gusmap = "~{sep=";" Gusmap_times_report}")

      #################################################################################
      # Table5 and 9: VariantEval  
      #################################################################################
      library(gsalib)
      df <- gsa.read.gatkreport("~{Freebayes_eval}")
      eval1 <- cbind(SNPCall = "Freebayes", seed = ~{seed}, depth = ~{depth}, df[["CompOverlap"]])
      count1 <- cbind(SNPCall = "Freebayes", seed = ~{seed}, depth = ~{depth}, df[["CountVariants"]])

      df <- gsa.read.gatkreport("~{GATK_eval}")
      eval2 <- cbind(SNPCall = "GATK", seed = ~{seed}, depth = ~{depth}, df[["CompOverlap"]])
      count2 <- cbind(SNPCall = "GATK", seed = ~{seed}, depth = ~{depth}, df[["CountVariants"]])
 
      df <- rbind(eval1, eval2)
      vroom_write(df, "data5_SNPCall_efficiency.tsv.gz", num_threads = ~{max_cores})

      df <- rbind(count1, count2)
      vroom_write(df, "data10_CountVariants.tsv.gz", num_threads = ~{max_cores})

      ##################################################################################
      # Table6: list of RDatas with name CountsFrom; seed; depth; SNPCall; GenoCall
      ##################################################################################

      default    <- str_split("~{sep=";" default_RDatas}", ";", simplify = T)
      SNPCaller  <- str_split("~{sep=";" SNPCaller_RDatas}", ";", simplify = T)
      Updog      <- str_split("~{sep=";" Updog_RDatas}", ";", simplify = T)
      Polyrad    <- str_split("~{sep=";" Polyrad_RDatas}", ";", simplify = T)
      Supermassa <- str_split("~{sep=";" Supermassa_RDatas}", ";", simplify = T)
      Gusmap <- str_split("~{sep=";" Gusmap_RDatas}", ";", simplify = T)

      RDatas_names <- c(default, SNPCaller, Updog, Polyrad, Supermassa, Gusmap)

      all_RDatas <- list()
      for(i in 1:length(RDatas_names)){
        map_temp <- load(RDatas_names[i])
        all_RDatas[[i]] <- get(map_temp)
      }
      all_RDatas <- unlist(all_RDatas, recursive = F)

      # Outputs
      vroom_write(errors_report, "data1_depths_geno_prob.tsv.gz", num_threads = ~{max_cores})
      vroom_write(maps_report, "data2_maps.tsv.gz", num_threads = ~{max_cores})
      vroom_write(filters_report, "data3_filters.tsv.gz", num_threads = ~{max_cores})
      vroom_write(times_report, "data4_times.tsv.gz", num_threads = ~{max_cores})

      gusmap_RDatas <- all_RDatas[grep("gusmap", names(all_RDatas))]
      RDatas <- all_RDatas[-grep("gusmap", names(all_RDatas))]

      # Converting OneMap sequencig objects to list. LargeList do not accept other class
      # Also because of this gusmap is separated, because the developers worked with enviroments, not classes

      for(i in 1:length(RDatas)){
        class(RDatas[[i]]) <- "list"
      }

      saveList(RDatas, file = "data6_RDatas.llo", append=FALSE, compress=TRUE)

      # LargeList package limits to 16 letters each character, therefore, the entire names are stored in a separated vector
      data_names <- as.data.frame(names(all_RDatas))
      vroom_write(data_names, "names.tsv.gz")
      save(gusmap_RDatas, file = "gusmap_RDatas.RData")

     RSCRIPT

  >>>

  runtime {
    docker:"cristaniguti/reads2map:0.0.1"
    preemptible: 3
    memory: "3 GB"
    cpu: 4
  }

  output {
    File data1_depths_geno_prob = "data1_depths_geno_prob.tsv.gz"
    File data2_maps = "data2_maps.tsv.gz"
    File data3_filters = "data3_filters.tsv.gz"
    File data4_times   = "data4_times.tsv.gz"
    File data5_SNPCall_efficiency = "data5_SNPCall_efficiency.tsv.gz"
    File data6_RDatas  = "data6_RDatas.llo"
    File data7_gusmap  = "gusmap_RDatas.RData"
    File data8_names   = "names.tsv.gz"
    File data10_counts  = "data10_CountVariants.tsv.gz"
  }
}

# Genotype calling with updog, polyrad and supermassa
# Exclusive for biallelic markers, input VCF file are normalized with -m-any
task ReGenotyping{
  input {
    String GenotypeCall_program
    File vcf_file
    String cross
    String parent1
    String parent2
    Int max_cores
  }

  command <<<
     R --vanilla --no-save <<RSCRIPT
       library(onemap)
       library(onemapUTILS)

       method <- "~{GenotypeCall_program}"
       cross <- "~{cross}"

       if(cross == "F1"){
          cross <- "outcross"
          f1 = NULL
       } else if (cross == "F2"){
          cross <- "f2 intercross"
          f1 = "F1"
       }

       out_vcf <- "regeno.vcf.gz"

        if (method == "updog") {
            out_onemap_obj <- updog_genotype(vcf="~{vcf_file}",
                                            vcf.par="AD",
                                            out_vcf = out_vcf,
                                            parent1="~{parent1}",
                                            parent2="~{parent2}",
                                            f1 = f1,
                                            crosstype= cross,
                                            recovering=TRUE,
                                            mean_phred=20,
                                            cores="~{max_cores}",
                                            depths=NULL,
                                            global_error = NULL,
                                            use_genotypes_errors = FALSE,
                                            use_genotypes_probs = TRUE)
        } else if (method == "supermassa") {
            out_onemap_obj <- supermassa_genotype(vcf="~{vcf_file}",
                                                  vcf.par="AD",
                                                  out_vcf = out_vcf,
                                                  parent1="~{parent1}",
                                                  parent2="~{parent2}",
                                                  crosstype= cross,
                                                  f1 = f1,
                                                  recovering=TRUE,
                                                  mean_phred=20,
                                                  cores="~{max_cores}",
                                                  depths=NULL,
                                                  global_error = NULL,
                                                  use_genotypes_errors = FALSE,
                                                  use_genotypes_probs = TRUE)
        } else if (method == "polyrad") {
            out_onemap_obj <- polyRAD_genotype(vcf="~{vcf_file}",
                                              parent1="~{parent1}",
                                              parent2="~{parent2}",
                                              f1 = f1,
                                              out_vcf = out_vcf,
                                              crosstype= cross,
                                              global_error = NULL,
                                              use_genotypes_errors = FALSE,
                                              use_genotypes_probs = TRUE,
                                              vcf.par="AD")
        }

        system("gunzip regeno.vcf.gz")

     RSCRIPT
  >>>

  runtime{
    docker:"cristaniguti/reads2map:0.0.1"
    # time:"20:00:00"
    # mem:"50GB"
    # cpu:20
    job_name: "ReGenotyping"
    node:"--nodes=1"
    mem:"--mem=30G"
    tasks:"--ntasks-per-node=15"
    time:"10:00:00"
  }

  output {
    File regeno_vcf = "regeno.vcf"
  }
}

task SetProbs{
  input{
    File vcf_file
    String cross
    String parent1
    String parent2
    String multiallelics
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
         f1 = "F1"
      }

      vcf <- read.vcfR("~{vcf_file}")
      save(vcf, file = "vcfR.RData")

      if("~{multiallelics}") only_biallelic = FALSE else only_biallelic = TRUE
          
      onemap.obj <- onemap_read_vcfR(vcfR.object = vcf,
                                     cross= cross,
                                     parent1="~{parent1}",
                                     parent2="~{parent2}",
                                     f1 = f1, only_biallelic = only_biallelic)

      #if(any(grepl("freeBayes", vcf@meta))) par <- "GL" else par <- "PL"

      probs <- extract_depth(vcfR.object=vcf,
                               onemap.object=onemap.obj,
                               vcf.par="GQ",
                               parent1="~{parent1}",
                               parent2="~{parent2}",
                               f1 = f1,
                               recovering=FALSE)

      probs_onemap_obj <- create_probs(input.obj = onemap.obj, genotypes_errors=probs)
      globalerror_onemap_obj <- create_probs(input.obj = onemap.obj, global_error = 0.05)

      save(probs_onemap_obj, file="probs_onemap_obj.RData")
      save(globalerror_onemap_obj, file="globalerror_onemap_obj.RData")

    RSCRIPT

  >>>
  runtime{
    docker:"cristaniguti/reads2map:0.0.1"
    # time:"10:00:00"
    # mem:"30GB"
    # cpu:1
    job_name: "SetProbs"
    node:"--nodes=1"
    mem:"--mem=20G"
    cpu:"--ntasks=1"
    time:"10:00:00"
  }

  output{
    File probs_onemap_obj = "probs_onemap_obj.RData"
    File globalerror_onemap_obj = "globalerror_onemap_obj.RData"
    File vcfR_obj = "vcfR.RData"
  }
}


task SetProbsDefault{
  input{
    File vcf_file
    String cross
    String parent1
    String parent2
    String multiallelics
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
        f1 = "F1"
      }

      vcf <- read.vcfR("~{vcf_file}")
      save(vcf, file = "vcfR.RData")
      
      if("~{multiallelics}") only_biallelic = FALSE else only_biallelic = TRUE

      onemap.obj <- onemap_read_vcfR(vcfR.object = vcf,
                                     cross= cross,
                                     parent1="~{parent1}",
                                     parent2="~{parent2}",
                                     f1 = f1, only_biallelic = only_biallelic)

      # if(any(grepl("freeBayes", vcf@meta))) par <- "GL" else par <- "PL"

      probs <- extract_depth(vcfR.object=vcf,
                               onemap.object=onemap.obj,
                               vcf.par="GQ",
                               parent1="~{parent1}",
                               parent2="~{parent2}",
                               f1 = f1,
                               recovering=FALSE)

      probs_onemap_obj <- create_probs(input.obj = onemap.obj, genotypes_errors=probs)
      globalerror_onemap_obj <- create_probs(input.obj = onemap.obj, global_error = 0.05)

      default_onemap_obj <- create_probs(input.obj = onemap.obj, global_error = 10^(-5))
      save(default_onemap_obj, file="default_onemap_obj.RData")
      save(probs_onemap_obj, file="probs_onemap_obj.RData")
      save(globalerror_onemap_obj, file="globalerror_onemap_obj.RData")

    RSCRIPT

  >>>
  runtime{
    docker:"cristaniguti/reads2map:0.0.1"
    # time:"10:00:00"
    # mem:"30GB"
    # cpu:1
    job_name: "SetProbsDefault"
    node:"--nodes=1"
    mem:"--mem=20G"
    cpu:"--ntasks=1"
    time:"10:00:00"
  }

  output{
    File probs_onemap_obj = "probs_onemap_obj.RData"
    File globalerror_onemap_obj = "globalerror_onemap_obj.RData"
    File default_onemap_obj = "default_onemap_obj.RData"
    File vcfR_obj = "vcfR.RData"
  }
}