version 1.0

import "./tasks/reads_simu.wdl" as sub

workflow SimulatedReads {

  input {
    Reference references
    Family family
    Sequencing sequencing
    Int number_of_families
    Int global_seed
  }

  # ProduceFamiliesSeeds just generates random seeds. It returns an
  # array of integers
  call ProduceFamiliesSeeds {
    input:
      global_seed= global_seed,
      number_of_families=number_of_families
  }

  # Here we generate Family objects on the fly, based on the values
  # from the family and the random seed of the previous task.
  scatter (seed in ProduceFamiliesSeeds.seeds) {
    Family fam =  {
      "cmBymb": family.cmBymb,
      "popsize": family.popsize,
      "enzyme1": sequencing.enzyme1,
      "enzyme2": sequencing.enzyme2,
      "seed": seed,
      "depth": sequencing.depth,
      "doses": family.doses,
      "ploidy": family.ploidy,
      "cross": family.cross,
      "multiallelics": sequencing.multiallelics
    }

    # Calling reads_simu for each seed
    call sub.reads_simu as ReadSimulations {
      input:
        references=references,
        family=fam,
        sequencing = sequencing
    }
  }

  call JointTables{
    input:
      data1=ReadSimulations.data1_depths_geno_prob,
      data2=ReadSimulations.data2_maps,
      data3=ReadSimulations.data3_filters,
      data5=ReadSimulations.data5_SNPcall_efficiency,
      data4=ReadSimulations.data4_times,
      depth=sequencing.depth,
      data6=ReadSimulations.data6_RDatas,
      data7=ReadSimulations.data7_gusmap,
      data8=ReadSimulations.data8_names,
      data9=ReadSimulations.simu_haplo,
      data10=ReadSimulations.multi_names
  }

  # Here you can reference outputs from the sub workflow. Remember that
  # it will be an array of the same type of the original.
  output {
    File results = JointTables.results
  }
}

task ProduceFamiliesSeeds {
  input {
    Int number_of_families
    Int global_seed
  }

  command <<<
    python <<CODE
    import random
    random.seed(~{global_seed})
    for x in range(~{number_of_families}):
        print(random.randint(1,101+x))
    CODE
  >>>

  runtime {
    docker: "python:3.7"
    time:"0:50:00"
    cpu:1
    mem:"--mem-per-cpu=14042"
  }

  output {
    Array[Int] seeds = read_lines(stdout())
  }
}


task JointTables{
  input {
    Array[File] data1
    Array[File] data2
    Array[File] data3
    Array[File] data4
    Array[File] data5
    Array[File] data6
    Array[File] data7
    Array[File] data8
    Array[File] data9
    Array[File] data10
    Int depth
  }

  command <<<

    R --vanilla --no-save <<RSCRIPT
    library(tidyverse)
    library(largeList)
    source("/opt/scripts/functions_simu.R")
    datas <- list()

    datas[[1]] <- c("~{sep=";" data1}")
    datas[[2]] <- c("~{sep=";" data2}")
    datas[[3]] <- c("~{sep=";" data3}")
    datas[[4]] <- c("~{sep=";" data4}")
    datas[[5]] <- c("~{sep=";" data5}")
    datas[[6]] <- c("~{sep=";" data6}")
    datas[[7]] <- c("~{sep=";" data7}")
    datas[[8]] <- c("~{sep=";" data8}")
    datas[[9]] <- c("~{sep=";" data9}")
    datas[[10]] <- c("~{sep=";" data10}")

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
      } else if(j == 10){
        multi_names_depth <- list()
        for(i in 1:length(datas[[j]])){
          multi_temp2 <- load(datas[[j]][i])
          multi_temp3 <- get(multi_temp2)
          multi_names_depth <- c(multi_names_depth, multi_temp3)
        }
      } else {
        for(i in 1:length(datas[[j]])){
          data_lst[[i]] <- readRDS(datas[[j]][i])
        }
        if(j == 8){
          dat <- do.call(c, data_lst)
        } else   dat <- do.call(rbind, data_lst)
        datas_up[[j]] <- dat
      }
    }

    result_list <- adapt2app(datas_up)

    saveRDS(result_list[[1]], file="data1.rds")
    saveRDS(result_list[[2]], file="data2.rds")
    saveRDS(result_list[[3]], file="data3.rds")
    saveRDS(result_list[[4]], file="data4.rds")
    saveRDS(result_list[[5]], file="data5.rds")
    saveRDS(datas_up[[9]], file="simu_haplo.rds")

    choices <- result_list[[6]]
    save(choices, file = "choices.RData")
    saveRDS(datas_up[[8]], file = "names.rds")
    save(multi_names_depth, file = "multi_names.RData")

    system("mkdir SimulatedReads_results_depth~{depth}")
    system("mv multi_names.RData gusmap_RDatas.RData sequences.llo data1.rds data2.rds data3.rds data4.rds data5.rds simu_haplo.rds choices.RData names.rds SimulatedReads_results_depth~{depth}")
    system("tar -czvf SimulatedReads_results_depth~{depth}.tar.gz SimulatedReads_results_depth~{depth}")

    RSCRIPT
  >>>

  runtime{
      docker:"cristaniguti/onemap_workflows"
      time:"03:00:00"
      cpu:1
      mem:"--mem-per-cpu=24042"
  }

  output{
    File results = "SimulatedReads_results_depth~{depth}.tar.gz"
  }
}
