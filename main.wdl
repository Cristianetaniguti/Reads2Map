version 1.0

import "./structs/reads_simuS.wdl"
import "./tasks/reads_simu.wdl" as sub

workflow main{

    input{
        ReferenceFasta references
        FamilyTemplate family_template
        Int number_of_families
    }

    # ProduceFamiliesSeeds just generates random seeds. It returns an
    # array of integers
    call ProduceFamiliesSeeds {
        input:
            number_of_families=number_of_families
    }

    # Here we generate Family objects on the fly, based on the values
    # from the family_template and the random seed of the previous task.
    scatter(seed in ProduceFamiliesSeeds.seeds) {
        Family fam =  {
            "cmBymb": family_template.cmBymb,
            "popsize": family_template.popsize,
            "enzyme": family_template.enzyme,
            "seed": seed,
            "depth": family_template.depth,
            "doses": family_template.doses,
            "ploidy": family_template.ploidy,
            "cross": family_template.cross
        }

        # Calling reads_simu for each seed
        call sub.reads_simu as ReadSimulations{
            input:
                references=references,
                family=fam
        }
    }

    call JointTables{
        input:
            data1  = ReadSimulations.data1_depths_geno_prob,
            data2  = ReadSimulations.data2_maps,
            data3  = ReadSimulations.data3_coverage,
            data4  = ReadSimulations.data4_filters,
            data5  = ReadSimulations.data5_SNPcall_efficiency,
    }

    # Here you can reference outputs from the sub workflow. Remember that
    # it will be an array of the same type of the original.
    output {
        File data1_depths_geno_prob = JointTables.data1_depths_geno_prob
        File data2_maps = JointTables.data2_maps
        File data3_coverage = JointTables.data3_coverage
        File data4_filters = JointTables.data4_filters
        File data5_SNPcall_efficiency = JointTables.data5_SNPcall_efficiency
    }
}

task ProduceFamiliesSeeds {
    input {
        Int number_of_families
    }

    command <<<
        python <<CODE
        import random
        for x in range(~{number_of_families}):
            print(random.randint(1,101))
        CODE
    >>>

    runtime {
        docker: "python:3.7"
    }

    output {
        Array[Int] seeds = read_lines(stdout())
    }
}



task JointTables{
    input{
        Array[File] data1 
        Array[File] data2
        Array[File] data3 
        Array[File] data4 
        Array[File] data5 
    }

    command <<<

        R --vanilla --no-save <<RSCRIPT

          datas <- list()

          datas[[1]] <- c("~{sep=";" data1}")
          datas[[2]] <- c("~{sep=";" data2}")
          datas[[3]] <- c("~{sep=";" data3}")
          datas[[4]] <- c("~{sep=";" data4}")
          datas[[5]] <- c("~{sep=";" data5}")

          datas <- lapply(datas, function(x) unlist(strsplit(x, ";")))
 
          data_lst <- list()
          for(j in 1:length(datas)){
              for(i in 1:length(datas[[j]])){
                  data_lst[[i]] <- readRDS(datas[[j]][i])
              }
          
            dat <- do.call(rbind, data_lst)
            saveRDS(dat, paste0("data",j,".rds"))
          }

        RSCRIPT
    >>>

    runtime{
        docker:"taniguti/onemap"
    }

    output{
       File data1_depths_geno_prob = "data1.rds"
       File data2_maps = "data2.rds"
       File data3_coverage = "data3.rds"
       File data4_filters = "data4.rds"
       File data5_SNPcall_efficiency = "data5.rds"
    }
}
