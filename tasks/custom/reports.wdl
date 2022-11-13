version 1.0


task JointReports{
  input{
    Array[File] SNPCaller
    Array[File] updog
    Array[File] polyrad
    Array[File] supermassa
    Array[File] gusmap
    Int max_cores
  }

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
    docker:"cristaniguti/reads2map:0.0.1"
    # time:"10:00:00"
    # mem:"80GB"
    # cpu:1
    job_name: "JointReports"
    node:"--nodes=1"
    mem:"--mem=20G"
    cpu:"--ntasks=4"
    time:"01:00:00"
  }

  output{
    File EmpiricalReads_results = "EmpiricalReads_results.tar.gz"
  }
}
