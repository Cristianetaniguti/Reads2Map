version 1.0

task vcf2onemap{
   input {
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
    runtime {
      docker:"cristaniguti/onemap_workflows"
      preemptible: 3
      memory: "2 GB"
      cpu:1
    }

    output{
      File onemap_obj = "~{SNPCall_program}_vcf_onemap.obj.RData"
      File vcfR_obj = "vcfR_obj.RData"
    }
}


task MultiVcf2onemap{
   input {
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
    runtime {
      docker:"cristaniguti/onemap_workflows"
      preemptible: 3
      memory: "3 GB"
      cpu: 1
    }

    output {
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

  runtime {
    docker: "cristaniguti/onemap_workflows"
    preemptible: 3
    memory:"3 GB"
    cpu:1
  }

  output {
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

  runtime {
    docker: "cristaniguti/onemap_workflows"
    time:"10:00:00"
    mem:"30GB"
    cpu:1
  }

  output {
    File filters_report = "filters_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
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

  runtime {
    docker: "cristaniguti/onemap_workflows"
    preemptible: 3
    memory: "8 GB"
    cpu: 4
  }

  output {
    File maps_report = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
    File maps_RData = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.RData"
    File times = "times_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
  }
}

task ErrorsReport {
  input {
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
    docker: "cristaniguti/onemap_workflows"
    preemptible: 3
    memory: "3 GB"
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
    preemptible: 3
    memory: "3 GB"
    cpu: 1
  }

  output {
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
       bam_vcf <- make_vcf(vcf_file, depths, allele_file, "~{program}_bam_vcf.vcf")

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
    mem:"30GB"
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
    mem:"30GB"
    cpu:4
  }

  output{
    File maps_report = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
    File maps_RData = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.RData"
    File times = "times_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
  }
}


task AddMultiallelics {
  input {
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

  runtime {
    docker:"cristaniguti/onemap_workflows"
    preemptible: 3
    memory: "8 GB"
    cpu: 4
  }

  output{
      File onemap_obj_both = "onemap_obj_both.RData"
  }
}


task JointReports{
  input {
    Array[File] default_RDatas
    Array[File] default_maps_report
    Array[File] default_filters_report
    Array[File] default_errors_report
    Array[File] default_times
    Array[File] multi_names
    Array[File] SNPCaller_RDatas
    Array[File] SNPCaller_maps_report
    Array[File] SNPCaller_filters_report
    Array[File] SNPCaller_errors_report
    Array[File] SNPCaller_times
    Array[File] Updog_RDatas
    Array[File] Updog_maps_report
    Array[File] Updog_filters_report
    Array[File] Updog_errors_report
    Array[File] Updog_times
    Array[File] Polyrad_RDatas
    Array[File] Polyrad_maps_report
    Array[File] Polyrad_filters_report
    Array[File] Polyrad_errors_report
    Array[File] Polyrad_times
    Array[File] Supermassa_RDatas
    Array[File] Supermassa_maps_report
    Array[File] Supermassa_filters_report
    Array[File] Supermassa_errors_report
    Array[File] Supermassa_times
    Array[File] Gusmap_RDatas
    Array[File] Gusmap_maps_report
    Array[File] Gusmap_times
    Int depth
    Int seed
    File gatk_ref_depth
    File gatk_ref_depth_bam
    File gatk_alt_depth
    File gatk_alt_depth_bam
    File freebayes_ref_depth_bam
    File freebayes_alt_depth_bam
    File freebayes_ref_depth
    File freebayes_alt_depth
  }

  command <<<
     R --vanilla --no-save <<RSCRIPT

      # Joint reports
      # I needed to split in groups because of R character limit size
      system("cat ~{sep= ' ' default_maps_report} ~{sep= ' ' SNPCaller_maps_report} > temp_map1")
      system("cat ~{sep= ' ' Updog_maps_report} ~{sep= ' ' Polyrad_maps_report} > temp_map2")
      system("cat ~{sep= ' ' Supermassa_maps_report}  ~{sep= ' ' Gusmap_maps_report} > temp_map3")
      system("cat temp_map1 temp_map2  temp_map3 > all_maps.txt")

      system("cat ~{sep= ' ' default_filters_report} ~{sep= ' ' SNPCaller_filters_report} > temp_filters1")
      system("cat ~{sep= ' ' Updog_filters_report} ~{sep= ' ' Polyrad_filters_report} > temp_filters2")
      system("cat ~{sep= ' ' Supermassa_filters_report} > temp_filters3")
      system("cat temp_filters1 temp_filters2 temp_filters3 > all_filters.txt")

      system("cat ~{sep= ' '  default_errors_report} ~{sep= ' ' SNPCaller_errors_report} > temp_errors1")
      system("cat ~{sep= ' ' Updog_errors_report} ~{sep= ' ' Polyrad_errors_report} > temp_errors2")
      system("cat ~{sep= ' ' Supermassa_errors_report} > temp_errors3")
      system("cat temp_errors1 temp_errors2 temp_errors3 > all_errors.txt")

      system("cat ~{sep= ' '  default_times} ~{sep= ' ' SNPCaller_times} > temp_times1")
      system("cat ~{sep= ' ' Updog_times} ~{sep= ' ' Polyrad_times} > temp_times2")
      system("cat ~{sep= ' ' Supermassa_times}  ~{sep= ' ' Gusmap_times} > temp_times3")
      system("cat temp_times1 temp_times2 temp_times3 > all_times.txt")

      # RDatas need to be load
      system("cp ~{sep= ' ' default_RDatas} ~{sep= ' ' SNPCaller_RDatas}  .")
      system("cp ~{sep= ' ' Updog_RDatas}  ~{sep= ' ' Polyrad_RDatas} .")
      system("cp ~{sep= ' ' Supermassa_RDatas} ~{sep= ' ' Gusmap_RDatas} .")

      # multiallelics names


      multi_temp <- c("~{sep=";" multi_names}")

      multi_temp <- unlist(strsplit(multi_temp, ";"))

      multi_names_seed <- list()
      for(i in 1:length(multi_temp)){
        multi_temp2 <- load(multi_temp[i])
        multi_temp3 <- get(multi_temp2)
        multi_names_seed <- c(multi_names_seed, multi_temp3)
      }

      save(multi_names_seed, file="multi_names.RData")

      library(tidyr)
      library(reshape2)
      library(largeList)

      depth <- ~{depth}
      seed <- ~{seed}

      all_errors <- read.table("all_errors.txt")

      # Adding allele counts to the report
      gatk.ref.depth <- read.table("~{gatk_ref_depth}")

      chr <- unique(sapply(strsplit(rownames(gatk.ref.depth), "_"), "[",1))

      all_errors[,5] <- paste0(chr, "_",all_errors[,5])
      colnames(all_errors) <- c("SNPCall", "GenoCall", "CountsFrom", "ind", "mks", "gabGT", "methGT", "A", "AB", "BA", "B")

      files <- list(c("~{gatk_ref_depth}", "~{gatk_alt_depth}"), c("~{gatk_ref_depth_bam}", "~{gatk_alt_depth_bam}"),
                    c("~{freebayes_ref_depth}", "~{freebayes_alt_depth}"), c("~{freebayes_ref_depth_bam}",
                                                                               "~{freebayes_alt_depth_bam}"))
      ref_alt <- c("ref", "alt")
      SNPCall <- rep(c("gatk", "freebayes"), each=4)
      CountsFrom <- rep(c("vcf", "bam"),2, each=2)

      ########################################################################################
      # Table1: GenoCall; mks; ind; SNPcall; CountsFrom; alt; ref; gabGT; methGT; A; AB; BA; B
      ########################################################################################

      # Binding depth information for each genotype
      z <- 1
      data1 <- data.frame()
      for(j in 1:length(files)){
        for(i in 1:2){
          if(i == 1) df_meth <- all_errors
          alleles <- read.table(files[[j]][i])
          alleles <- cbind(mks=rownames(alleles), alleles)
          alleles <- melt(alleles, id.vars = c("mks"))
          colnames(alleles) <- c("mks", "ind", ref_alt[i])
          alleles <- cbind(SNPCall=SNPCall[z], CountsFrom = CountsFrom[z], alleles)
          z <- z+ 1
          df_meth <- merge(df_meth, alleles, by = c("SNPCall", "CountsFrom", "ind", "mks"))

        }
        data1 <- rbind(data1, df_meth)
      }

      # Adding seed and depth info
      all_errors <- cbind(seed, depth, data1)

      ########################################################
      # Table2: seed; CountsFrom; ErrorProb; SNPcall; MK; rf; phases; simulated_phases
      ########################################################
      # Add seed and mean depth information
      all_maps <- read.table("all_maps.txt")
      all_maps <- cbind(seed=seed, depth=depth, all_maps)

      ##########################################################################
      # Table3: CountsFrom; seed; SNPcall; GenoCall; n_mks; distorted; redundant
      ##########################################################################
      all_filters <- read.table("all_filters.txt")
      all_filters <- cbind(seed=seed, depth=depth, all_filters)

      ###########################################################################
      # Table4: CountsFrom; seed; SNPcall; GenoCall
      ###########################################################################
      all_times <- read.table("all_times.txt", stringsAsFactors = F)

      temp <- strsplit(all_times[,2], "_")
      temp <- do.call(rbind, temp)
      temp <- temp[,-1]
      colnames(temp) <- c("SNPCall", "CountsFrom", "GenoCall", "fake")
      all_times <- cbind(seed=seed, depth=depth, temp, time= all_times[,3])

      ###########################################################################
      # Table6: list of RDatas with name CountsFrom; seed; depth; SNPcall; GenoCall
      ###########################################################################
      Genocall <- c("default", "SNPCaller", "updog", "supermassa", "polyrad", "gusmap",
                    "default0.05", "updog0.05", "supermassa0.05", "polyrad0.05")
      fake <- c(TRUE, FALSE)

      CountsFrom <- c("vcf", "bam")
      SNPCall <- c("gatk", "freebayes")
      df <- data.frame(SNPCall, CountsFrom, Genocall)
      df <- tidyr::expand(df, SNPCall, CountsFrom, Genocall)
      df <- as.data.frame(df)
      df <- df[-which((df[,3] == "default" | df[,3] == "default0.05" | df[,3] == "SNPCaller" ) & df[,2] == "bam"),]
      RDatas_names <- paste0("map_",df[,1],"_",df[,2], "_", df[,3],".RData")

      all_RDatas <- list()
      for(i in 1:length(RDatas_names)){
         map_temp <- load(RDatas_names[i])
         all_RDatas[[i]] <- get(map_temp)
      }
      all_RDatas <- unlist(all_RDatas, recursive = F)

      all_names <- names(all_RDatas)
      new_names <- paste0(seed, "_", depth, "_", all_names)
      names(all_RDatas) <- new_names

      # Outputs
      saveRDS(all_errors, file = "data1_depths_geno_prob.rds")
      saveRDS(all_maps, file = "data2_maps.rds")
      saveRDS(all_filters, file = "data3_filters.rds")
      saveRDS(all_times, file= "data4_times.rds")

      gusmap_RDatas <- all_RDatas[grep("gusmap", names(all_RDatas))]
      RDatas <- all_RDatas[-grep("gusmap", names(all_RDatas))]

      # Converting OneMap sequencig objects to list. LargeList do not accept other class
      # Also because of this gusmap is separated, because the developers worked with enviroments, not classes

      for(i in 1:length(RDatas)){
        class(RDatas[[i]]) <- "list"
      }

      saveList(RDatas, file = "data6_RDatas.llo", append=FALSE, compress=TRUE)

      # LargeList package limits to 16 letters each character, therefore, the entire names are stored in a separated vector
      saveRDS(new_names, file = "names.rds")
      save(gusmap_RDatas, file = "gusmap_RDatas.RData")

     RSCRIPT

  >>>

  runtime {
    docker:"cristaniguti/onemap_workflows"
    preemptible: 3
    memory: "4 GB"
    cpu: 1
  }

  output {
    File data1_depths_geno_prob = "data1_depths_geno_prob.rds"
    File data2_maps = "data2_maps.rds"
    File data3_filters = "data3_filters.rds"
    File data4_times = "data4_times.rds"
    File data6_RDatas = "data6_RDatas.llo"
    File data7_gusmap = "gusmap_RDatas.RData"
    File data8_names = "names.rds"
    File multi_names2 = "multi_names.RData"
  }
}
