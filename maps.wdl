version 1.0

import "./structs/maps.wdl"

workflow maps {
  input {
    outputF2 outputf2
  }
  
  Array[String] methodNames = read_lines(outputf2.methods)  
  
  scatter (methodName in methodNames){
    call default{
      input:
        vcfRs = outputf2.vcfRs,
        tot_mks = outputf2.tot_mks,
        methodName = methodName
    }
  }
}

task default {
  input {
    File vcfRs
    File tot_mks
    String methodName
  }
  output{
    File filters = "~{methodName}_filters.txt"
    File map_df = "~{methodName}_map_df.txt"
  }
  command <<<

        R --vanilla --no-save <<RSCRIPT
        load("~{vcfRs}")
        library(onemap)
        tot_mks <- read.table("~{tot_mks}")

        if("~{methodName}" == "stacks"){
        df <- onemap_read_vcfR(vcfR.object=vcfRs[["~{methodName}"]], 
                                cross="f2 intercross", 
                                parent1="P1_rg", 
                                parent2="P2_rg", 
                                f1="F1_rg")
        } else {
          df <- onemap_read_vcfR(vcfR.object=vcfRs[["~{methodName}"]], 
                                cross="f2 intercross", 
                                parent1="P1", 
                                parent2="P2", 
                                f1="F1")
        }
        ## Filters        
        n.mk <- df[[3]]
        segr <- test_segregation(df)
        distorted <- length(select_segreg(segr, distorted = T))
        bins <- find_bins(df)
        nbins <- length(bins[[1]])

        filters <- data.frame("n_markers"= n.mk, 
                              "distorted_markers"=distorted, 
                              "redundant_markers"=n.mk-nbins)
        write.table(filters, file = paste0("~{methodName}","_filters.txt"), row.names = F, quote = F)

        ## Map
        twopts <- rf_2pts(df)
        true.mks <- which(df[[9]] %in% tot_mks[,2])
        seq.true <- make_seq(twopts, true.mks)
        map.df <- map(seq.true)
        map.info <- data.frame("mk.name"= colnames(df[[1]])[map.df[[1]]], 
                              "pos" = df[[9]][map.df[[1]]],
                              "rf" = c(0,cumsum(haldane(map.df[[3]]))))
        write.table(map.info, file = paste0("~{methodName}", "_map_df.txt"), row.names = F, quote = F)


        RSCRIPT
  >>>

  runtime{
    docker:"onemap:v1"
  }
}
