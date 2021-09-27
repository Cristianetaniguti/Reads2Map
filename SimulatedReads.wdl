version 1.0

import "./tasks/SimulatedSingleFamily.wdl" as sub

workflow SimulatedReads {

  input {
    Reference references
    Family family
    Sequencing sequencing
    Int number_of_families
    Int global_seed
    Int max_cores
    String? filters
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
    call sub.SimulatedSingleFamily {
      input:
        references=references,
        family=fam,
        sequencing = sequencing,
        max_cores = max_cores,
        filters = filters
    }
  }

  call JointTables {
    input:
      data1_depths_geno_prob   = SimulatedSingleFamily.data1_depths_geno_prob,
      data2_maps               = SimulatedSingleFamily.data2_maps,
      data3_filters            = SimulatedSingleFamily.data3_filters,
      data5_SNPCall_efficiency = SimulatedSingleFamily.data5_SNPCall_efficiency,
      data4_times              = SimulatedSingleFamily.data4_times,
      data6_RDatas             = SimulatedSingleFamily.data6_RDatas,
      data7_gusmap             = SimulatedSingleFamily.data7_gusmap,
      data8_names              = SimulatedSingleFamily.data8_names,
      data9_simu_haplo         = SimulatedSingleFamily.simu_haplo,
      data10_counts            = SimulatedSingleFamily.data10_counts,
      depth                    = sequencing.depth,
      plots                    = SimulatedSingleFamily.Plots
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
    preemptible: 3
    cpu: 1
    memory: "1 GB"
  }

  output {
    Array[Int] seeds = read_lines(stdout())
  }
}


task JointTables{
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
    Array[File] plots
    Int depth        
  }

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
            simu_haplo.tsv.gz  names.tsv.gz ~{sep=" " plots} SimulatedReads_results_depth~{depth}")
    system("tar -czvf SimulatedReads_results_depth~{depth}.tar.gz SimulatedReads_results_depth~{depth}")

    RSCRIPT
  >>>

  runtime {
      docker:"cristaniguti/reads2map:0.0.1"
      preemptible: 3
      cpu: 1
      memory: "3 GB"
  }

  output {
    File results = "SimulatedReads_results_depth~{depth}.tar.gz"
  }
}
