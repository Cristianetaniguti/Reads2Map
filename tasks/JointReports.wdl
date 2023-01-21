version 1.0

task JointReports{
  input{
    Array[File?] SNPCaller
    Array[File?] updog
    Array[File?] polyrad
    Array[File?] supermassa
    Array[File?] gusmap
    Int max_cores
  }

  Int disk_size = ceil(size(SNPCaller, "GiB") * 1.5 + size(updog, "GiB") * 1.5 + size(polyrad, "GiB") * 1.5 + size(supermassa, "GiB") * 1.5 + size(gusmap, "GiB") * 1.5)
  Int memory_size = 4000

  command <<<
     R --vanilla --no-save <<RSCRIPT

     # packages
      library(tidyr)
      library(stringr)
      library(vroom)
      library(largeList)

      SNPCaller  <- str_split("~{sep=";" SNPCaller}", ";", simplify = T)
      updog      <- str_split("~{sep=";" updog}", ";", simplify = T)
      polyrad    <- str_split("~{sep=";" polyrad}", ";", simplify = T)
      supermassa <- str_split("~{sep=";" supermassa}", ";", simplify = T)
      gusmap <- str_split("~{sep=";" gusmap}", ";", simplify = T)

      files <- list(SNPCaller, updog, polyrad, supermassa, gusmap)

      path_dir <- tempdir()
      system(paste0("mkdir ", paste0(path_dir, c("/maps", "/filters", "/errors", "/times", "/RDatas"), collapse = " ")))
      for(i in 1:length(files)){
        for(j in 1:length(files[[i]])){
          untar(files[[i]][[j]], exdir = path_dir)
          list_files <- untar(files[[i]][[j]], exdir = path_dir, list = T)
          system(paste0("mv ",path_dir, "/",list_files[1], "*_map_report.tsv.gz ", path_dir, "/maps"))
          system(paste0("mv ",path_dir, "/",list_files[1], "*_times_report.tsv.gz ", path_dir, "/times"))
          system(paste0("mv ",path_dir, "/",list_files[1], "*.RData ", path_dir, "/RDatas"))
          if(!grepl("gusmap", list_files[1])){
            system(paste0("mv ",path_dir, "/",list_files[1], "*_filters_report.tsv.gz ", path_dir, "/filters"))
            system(paste0("mv ",path_dir, "/",list_files[1], "*_errors_report.tsv.gz ", path_dir, "/errors"))
          }
        }
      }

      files <- system(paste0("ls ", path_dir, "/maps/"), intern = T)
      files <- paste0(path_dir, "/maps/", files)
      maps <- vroom(files, num_threads = ~{max_cores})

      files <- system(paste0("ls ", path_dir, "/filters/"), intern = T)
      files <- paste0(path_dir, "/filters/", files)
      filters <- vroom(files, num_threads = ~{max_cores})

      files <- system(paste0("ls ", path_dir, "/errors/"), intern = T)
      files <- paste0(path_dir, "/errors/", files)
      errors <- vroom(files, num_threads = ~{max_cores})

      files <- system(paste0("ls ", path_dir, "/times/"), intern = T)
      files <- paste0(path_dir, "/times/", files)
      times <- vroom(files, num_threads = ~{max_cores})

      files <- system(paste0("ls ", path_dir, "/RDatas/"), intern = T)
      files <- paste0(path_dir, "/RDatas/", files)

      all_RDatas <- list()
      for(i in 1:length(files)){
        map_temp <- load(files[i])
        all_RDatas[[i]] <- get(map_temp)
      }

      names(all_RDatas) <- basename(files)
      gusmap_RDatas <- all_RDatas[grep("gusmap", names(all_RDatas))]
      RDatas <- all_RDatas[-grep("gusmap", names(all_RDatas))]

      #   # Converting onemap sequencig objects to list. LargeList do not accept other class
      #   # Also because of this gusmap is separated, because the developers worked with enviroments, not classes

      for(i in 1:length(RDatas)){
        class(RDatas[[i]]) <- "list"
      }

      saveList(RDatas, file = "sequences_emp.llo", append=FALSE, compress=TRUE)

      new_names <- names(all_RDatas)
      vroom_write(as.data.frame(new_names), "names.tsv.gz")
      save(gusmap_RDatas, file = "gusmap_RDatas.RData")

      # Outputs
      vroom_write(errors, "data1_depths_geno_prob.tsv.gz", num_threads = ~{max_cores})
      vroom_write(maps, "data2_maps.tsv.gz", num_threads = ~{max_cores})
      vroom_write(filters, "data3_filters.tsv.gz", num_threads = ~{max_cores})
      vroom_write(times, "data4_times.tsv.gz", num_threads = ~{max_cores})

      system("mkdir EmpiricalReads_results")
      system("mv gusmap_RDatas.RData sequences_emp.llo data1_depths_geno_prob.tsv.gz data2_maps.tsv.gz data3_filters.tsv.gz data4_times.tsv.gz names.tsv.gz EmpiricalReads_results")
      system("tar -czvf EmpiricalReads_results.tar.gz EmpiricalReads_results")

     RSCRIPT

  >>>

  runtime{
    docker:"cristaniguti/reads2map:0.0.4"
    cpu: max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "JointReports"
    mem:"~{memory_size}M"
    time:"01:40:00"
  }

  output{
    File EmpiricalReads_results = "EmpiricalReads_results.tar.gz"
  }
}

task JointReportsSimu {
  input {
    Array[File] SNPCaller
    Array[File] updog
    Array[File] polyrad
    Array[File] supermassa
    Array[File] gusmap_files
    Array[File] multiallelics_file
    File Freebayes_eval
    File GATK_eval
    Int max_cores
    Int seed
    Int depth
  }

  Int disk_size = ceil(size(SNPCaller, "GiB") * 1.5 + size(updog, "GiB") * 1.5 + size(polyrad, "GiB") * 1.5 + size(supermassa, "GiB") * 1.5 + size(gusmap_files, "GiB") * 1.5 + size(multiallelics_file, "GiB") * 1.5)
  Int memory_size = 5000

  command <<<
     R --vanilla --no-save <<RSCRIPT

      # packages
      library(tidyr)
      library(stringr)
      library(vroom)
      library(largeList)
      library(vcfR)

      SNPCaller  <- str_split("~{sep=";" SNPCaller}", ";", simplify = T)
      updog      <- str_split("~{sep=";" updog}", ";", simplify = T)
      polyrad    <- str_split("~{sep=";" polyrad}", ";", simplify = T)
      supermassa <- str_split("~{sep=";" supermassa}", ";", simplify = T)
      gusmap <- str_split("~{sep=";" gusmap_files}", ";", simplify = T)

      files <- list(SNPCaller, updog, polyrad, supermassa, gusmap)

      direc <- c("/maps/", "/filters/", "/errors/", "/times/", "/RDatas/")

      path_dir <- tempdir()
      system(paste0("mkdir ", paste0(path_dir, direc, collapse = " ")))
      for(i in 1:length(files)){
        for(j in 1:length(files[[i]])){
            untar(files[[i]][[j]], exdir = path_dir)
            list_files <- untar(files[[i]][[j]], exdir = path_dir, list = T)
            system(paste0("mv ",path_dir, "/",list_files[1], "*_map_report.tsv.gz ", path_dir, "/maps"))
            system(paste0("mv ",path_dir, "/",list_files[1], "*_times_report.tsv.gz ", path_dir, "/times"))
            system(paste0("mv ",path_dir, "/",list_files[1], "*.RData ", path_dir, "/RDatas"))
            if(!grepl("gusmap", list_files[1])){
              system(paste0("mv ",path_dir, "/",list_files[1], "*_filters_report.tsv.gz ", path_dir, "/filters"))
              system(paste0("mv ",path_dir, "/",list_files[1], "*_errors_report.tsv.gz ", path_dir, "/errors"))
            }
        }
      }

      direc_tsv <- direc[-5]
      tsvs <- list()
      for(i in 1:length(direc_tsv)){
        files <- system(paste0("ls ", path_dir, direc_tsv[i]), intern = T)
        files <- paste0(path_dir, direc_tsv[i], files)
        tsvs[[i]] <- vroom(files, num_threads = ~{max_cores})
      }

      # Add multiallelics tag
      vcf_multi <- str_split("~{sep = ";" multiallelics_file}", ";", simplify = T)

      vcfs <- lapply(as.list(vcf_multi), read.vcfR)
      multi_names_seed <- lapply(vcfs, function(x) paste0(x@fix[,1], "_", x@fix[,2]))

      snpcall_names <- rep(NA, length(vcfs))
      snpcall_names[which(sapply(vcfs, function(x) any(grep("gatk",x@meta))))] <- "gatk"
      snpcall_names[which(sapply(vcfs, function(x) any(grep("freebayes",x@meta))))] <- "freebayes"

      for(i in 1:length(snpcall_names)){
        tsvs[[1]][,"real.mks"][which(tsvs[[1]][,"SNPCall"] == snpcall_names[i] &
                                      tsvs[[1]][,"mk.name"] %in% multi_names_seed[[i]])] <- "multiallelic"
      }

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

      rdatas_files <- paste0(path_dir, "/RDatas/",list.files(paste0(path_dir, "/RDatas/")))

      all_RDatas <- list()
      for(i in 1:length(rdatas_files)){
        map_temp <- load(rdatas_files[i])
        all_RDatas[[i]] <- get(map_temp)
      }
      all_RDatas <- unlist(all_RDatas, recursive = F)

      gusmap_RDatas <- all_RDatas[grep("gusmap", names(all_RDatas))]
      RDatas <- all_RDatas[-grep("gusmap", names(all_RDatas))]

      #   # Converting onemap sequencig objects to list. LargeList do not accept other class
      #   # Also because of this gusmap is separated, because the developers worked with enviroments, not classes

      for(i in 1:length(RDatas)){
        class(RDatas[[i]]) <- "list"
      }

      saveList(RDatas, file = "data6_RDatas.llo", append=FALSE, compress=TRUE)

      new_names <- names(all_RDatas)
      vroom_write(as.data.frame(new_names), "names.tsv.gz")
      save(gusmap_RDatas, file = "gusmap_RDatas.RData")

      # Outputs
      vroom_write(tsvs[[3]], "data1_depths_geno_prob.tsv.gz", num_threads = ~{max_cores})
      vroom_write(tsvs[[1]], "data2_maps.tsv.gz", num_threads = ~{max_cores})
      vroom_write(tsvs[[2]], "data3_filters.tsv.gz", num_threads = ~{max_cores})
      vroom_write(tsvs[[4]], "data4_times.tsv.gz", num_threads = ~{max_cores})

     RSCRIPT
  >>>

  runtime {
    docker:"cristaniguti/reads2map:0.0.4"
    cpu: max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "JointReports"
    mem:"~{memory_size}M"
    time:"01:40:00"
  }

  meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Merge reports resulted from evaluations with different pipelines."
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

task JointTablesSimu{
  input {
    Array[File] data1_depths_geno_prob
    Array[File] data2_maps
    Array[File] data3_filters
    Array[File] data5_SNPCall_efficiency
    Array[File] data4_times
    Array[File] data6_RDatas
    Array[File] data7_gusmap
    Array[File] data8_names
    Array[File] data9_simu_haplo
    Array[File] data10_counts
    Array[File?] plots
    Array[File] positions
    Int depth
  }

  Int disk_size = ceil(size(data1_depths_geno_prob, "GiB") * 2 + size(data2_maps, "GiB") + size(data3_filters, "GiB") + size(data5_SNPCall_efficiency, "GiB") + size(data4_times, "GiB") + size(data6_RDatas, "GiB") + size(data7_gusmap, "GiB") + size(data8_names, "GiB") + size(data9_simu_haplo, "GiB") + size(data10_counts, "GiB"))
  Int memory_size = 7000

  command <<<

    R --vanilla --no-save <<RSCRIPT
    library(tidyverse)
    library(largeList)
    library(vroom)

    datas <- list()

    datas[[1]] <- c("~{sep=";" data1_depths_geno_prob  }")
    datas[[2]] <- c("~{sep=";" data2_maps              }")
    datas[[3]] <- c("~{sep=";" data3_filters           }")
    datas[[4]] <- c("~{sep=";" data5_SNPCall_efficiency}")
    datas[[5]] <- c("~{sep=";" data4_times             }")
    datas[[6]] <- c("~{sep=";" data6_RDatas            }")
    datas[[7]] <- c("~{sep=";" data7_gusmap            }")
    datas[[8]] <- c("~{sep=";" data8_names             }")
    datas[[9]] <- c("~{sep=";" data9_simu_haplo        }")
    datas[[10]] <- c("~{sep=";" data10_counts        }")

    datas <- lapply(datas, function(x) unlist(strsplit(x, ";")))

    Rdata_lst <- data_lst <- datas_up <- list()
    for(j in 1:length(datas)){
      if(j == 6){
        for(i in 1:length(datas[[j]])){
          temp <- readList(datas[[j]][i])
          if(i == 1){
            saveList(temp, file="sequences.llo", append = F, compress = T)
          } else {
            saveList(temp, file="sequences.llo", append = T, compress = T)
          }
        }
      } else  if(j == 7){
        for(i in 1:length(datas[[j]])){
          temp <- load(datas[[j]][i])
          Rdata_lst[[i]] <- get(temp)
        }
        Rdatas <- do.call(c, Rdata_lst)
        save(Rdatas, file = "gusmap_RDatas.RData")
      } else {
        for(i in 1:length(datas[[j]])){
          data_lst[[i]] <- vroom(datas[[j]][i], delim = "\t")
        }
        if(j == 8){
          dat <- do.call(c, data_lst)
        } else   dat <- do.call(rbind, data_lst)
        datas_up[[j]] <- dat
      }
    }

    vroom_write(datas_up[[1]], "data1_depths_geno_prob.tsv.gz")
    vroom_write(datas_up[[2]], "data2_maps.tsv.gz")
    vroom_write(datas_up[[3]], "data3_filters.tsv.gz")
    vroom_write(datas_up[[5]], "data4_times.tsv.gz")
    vroom_write(datas_up[[4]], "data5_SNPCall_efficiency.tsv.gz")
    vroom_write(datas_up[[9]], "simu_haplo.tsv.gz")
    vroom_write(datas_up[[10]], "data10_counts.tsv.gz")

    data.names <- as.data.frame(datas_up[[8]])
    print(data.names)
    vroom_write(data.names, "names.tsv.gz")

    system("mkdir SimulatedReads_results_depth~{depth}")
    system("mv gusmap_RDatas.RData sequences.llo data1_depths_geno_prob.tsv.gz \
            data2_maps.tsv.gz data3_filters.tsv.gz data4_times.tsv.gz data5_SNPCall_efficiency.tsv.gz data10_counts.tsv.gz \
            simu_haplo.tsv.gz  names.tsv.gz ~{sep=" " plots} ~{sep=" " positions} SimulatedReads_results_depth~{depth}")
    system("tar -czvf SimulatedReads_results_depth~{depth}.tar.gz SimulatedReads_results_depth~{depth}")

    RSCRIPT
  >>>

  runtime {
    docker:"cristaniguti/reads2map:0.0.4"
    cpu: 1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "JointTables"
    mem:"~{memory_size}M"
    time:"01:40:00"     
  }

  output {
    File results = "SimulatedReads_results_depth~{depth}.tar.gz"
  }
}

task JointReportsPoly{
  input{
    Array[File?] SNPCaller
    Array[File?] updog
    Array[File?] polyrad
    Array[File?] supermassa
  }

  Int disk_size = ceil(size(SNPCaller, "GiB") * 1.5 + size(updog, "GiB") * 1.5 + size(polyrad, "GiB") * 1.5 + size(supermassa, "GiB") * 1.5)
  Int memory_size = 4000

  command <<<

     snpcaller=(~{sep=" " SNPCaller})
     updog=(~{sep=" " updog})
     polyrad=(~{sep=" " polyrad})
     supermassa=(~{sep=" " supermassa})

     mkdir results_all

     for index in ${!snpcaller[*]}; 
       tar -xvf ${snpcaller[$index]}
       mv results/* results_all
       rm -r results
     done

     for index in ${!updog[*]}; 
       tar -xvf ${updog[$index]}
       mv results/* results_all
       rm -r results     
     done

     for index in ${!polyrad[*]}; 
       tar -xvf ${polyrad[$index]}
       mv results/* results_all
       rm -r results     
     done

     for index in ${!supermassa[*]}; 
       tar -xvf ${supermassa[$index]}
       mv results/* results_all
       rm -r results
     done

    tar -czvf EmpiricalReads_results_poly.tar.gz results_all

  >>>

  runtime{
    docker:"ubuntu:20.04"
    cpu: 1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "JointReports"
    mem:"~{memory_size}M"
    time:"01:40:00"
  }

  output{
    File EmpiricalReads_results = "EmpiricalReads_results_poly.tar.gz"
  }
}