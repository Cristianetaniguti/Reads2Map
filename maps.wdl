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
        tot_mks = outputf2.tot_mks,
        methodName = methodName,
        simu_vcf = outputf2.simu_vcf,
        gatkVCF = outputf2.gatkVCF,
        freebayesVCF = outputf2.freebayesVCF,
        stacksVCF = outputf2.stacksVCF
    }
  }
}

task default {
  input {
    File tot_mks
    String methodName
    File simu_vcf
    File gatkVCF
    File freebayesVCF
    File stacksVCF
  }
  output{
    File filters = "~{methodName}_filters.txt"
    File map_df = "~{methodName}_map_df.txt"
    File map_GQ = "~{methodName}_map_GQ.txt"
    File GQ_error_info = "~{methodName}_error_info.txt"
    File updog_error_info = "~{methodName}_updog_error_info.txt"
  }
  command <<<

        R --vanilla --no-save <<RSCRIPT
        library(onemap)
        library(updog)
        library(reshape2)
        library(vcfR)

        tot_mks <- read.table("~{tot_mks}")

        if("~{methodName}" == "gatk"){
          vcf <- read.vcfR("~{gatkVCF}")
        }
        if("~{methodName}" == "freebayes"){
          vcf <- read.vcfR("~{freebayesVCF}")
        }
        if("~{methodName}" == "stacks"){
          vcf <- read.vcfR("~{stacksVCF}")
        }

        if("~{methodName}" == "stacks"){
        df <- onemap_read_vcfR(vcfR.object=vcf, 
                                cross="f2 intercross", 
                                parent1="P1_rg", 
                                parent2="P2_rg", 
                                f1="F1_rg")
        } else {
          df <- onemap_read_vcfR(vcfR.object=vcf, 
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

        ## Maps default
        twopts <- rf_2pts(df)
        true.mks <- which(df[[9]] %in% tot_mks[,2])
        seq.true <- make_seq(twopts, true.mks)
        map.df <- map(seq.true)
        map.info <- data.frame("mk.name"= colnames(df[[1]])[map.df[[1]]], 
                              "pos" = df[[9]][map.df[[1]]],
                              "rf" = c(0,cumsum(haldane(map.df[[3]]))))
        write.table(map.info, file = paste0("~{methodName}", "_map_df.txt"), row.names = F, quote = F)

        ## extract GQ
        
        if("~{methodName}" == "stacks"){
          aval.gq <- extract_depth(vcfR.object= vcf,
                                onemap.object = df,
                                vcf.par = "GQ",
                                parent1 = "P1_rg",
                                parent2 = "P2_rg",
                                f1="F1_rg",
                                recovering = FALSE)
        } else {
          aval.gq <- extract_depth(vcfR.object= vcf,
                                onemap.object = df,
                                vcf.par = "GQ",
                                parent1 = "P1",
                                parent2 = "P2",
                                f1="F1",
                                recovering = FALSE)
        }

        ## Maps GQ
        twopts <- rf_2pts(aval.gq)
        true.mks <- which(aval.gq[[9]] %in% tot_mks[,2])
        seq.true <- make_seq(twopts, true.mks)
        map.df <- map(seq.true)
        map.info <- data.frame("mk.name"= colnames(aval.gq[[1]])[map.df[[1]]], 
                              "pos" = aval.gq[[9]][map.df[[1]]],
                              "rf" = c(0,cumsum(haldane(map.df[[3]]))))
        write.table(map.info, file = paste0("~{methodName}", "_map_GQ.txt"), row.names = F, quote = F)

        # Errors info tab
        simu <- read.vcfR("~{simu_vcf}")
        gab <- onemap_read_vcfR(vcfR.object = simu,
                                cross = "f2 intercross",
                                parent1 = "P1",
                                parent2 = "P2",
                                f1 = "F1")
        
        pos <- which(gab[[9]] %in% aval.gq[[9]])
        pos.inv <- which(aval.gq[[9]] %in% gab[[9]])
        
        gab.pos <- gab[[9]][pos]
        gab.geno <- gab[[1]][,pos]
        colnames(gab.geno) <- gab.pos
        gab.geno <-melt(gab.geno)
        colnames(gab.geno) <- c("MK", "POS", "gabGT")
        
        meth.geno <- aval.gq[[1]][,pos.inv]
        meth.error <- aval.gq[[10]][,pos.inv]
        meth.pos <- aval.gq[[9]][pos.inv]
        colnames(meth.error) <- colnames(meth.geno) <- meth.pos
        meth.geno <- melt(meth.geno)
        colnames(meth.geno) <- c("MK", "POS", "methGT")
        
        meth.error <- melt(meth.error)
        colnames(meth.error) <- c("MK", "POS", "methError")
        
        error.info <- merge(gab.geno, meth.geno)
        error.info <- merge(error.info, meth.error)
        error.info <- error.info[order(error.info[,2], as.character(error.info[,1])),]
        
         write.table(error.info, file = paste0("~{methodName}", "_error_info.txt"), row.names = F, quote = F)

        ## Updog
        
        if("~{methodName}" == "stacks"){
          updog.aval <- updog_error(vcfR.object=vcf,
                                    onemap.object = df,
                                    vcf.par = "AD",
                                    parent1 = "P1_rg",
                                    parent2 = "P2_rg",
                                    f1="F1_rg",
                                    recovering = TRUE,
                                    mean_phred = 20,
                                    cores = 6,
                                    depths = NULL)
        } else {
          updog.aval <- updog_error(vcfR.object=vcf,
                                    onemap.object = df,
                                    vcf.par = "AD",
                                    parent1 = "P1",
                                    parent2 = "P2",
                                    f1="F1",
                                    recovering = TRUE,
                                    mean_phred = 20,
                                    cores = 6,
                                    depths = NULL)
        }

        ## Filters        
        n.mk <- updog.aval[[3]]
        segr <- test_segregation(updog.aval)
        distorted <- length(select_segreg(segr, distorted = T))
        bins <- find_bins(updog.aval)
        nbins <- length(bins[[1]])

        filters <- data.frame("n_markers"= n.mk, 
                              "distorted_markers"=distorted, 
                              "redundant_markers"=n.mk-nbins)
        write.table(filters, file = paste0("~{methodName}","_filters_updog.txt"), row.names = F, quote = F)

        # Maps updog
        twopts <- rf_2pts(updog.aval)
        true.mks <- which(updog.aval[[9]] %in% tot_mks[,2])
        seq.true <- make_seq(twopts, true.mks)
        map.df <- map(seq.true)
        map.info <- data.frame("mk.name"= colnames(updog.aval[[1]])[map.df[[1]]], 
                               "pos" = updog.aval[[9]][map.df[[1]]],
                               "rf" = c(0,cumsum(haldane(map.df[[3]]))))
        write.table(map.info, file = paste0("~{methodName}", "_map_updog.txt"), row.names = F, quote = F)


        # Errors info
        pos.inv <- which(updog.aval[[9]] %in% gab[[9]])
        
        meth.geno <- updog.aval[[1]][,pos.inv]
        meth.error <- updog.aval[[10]][,pos.inv]
        meth.pos <- updog.aval[[9]][pos.inv]
        colnames(meth.error) <- colnames(meth.geno) <- meth.pos
        meth.geno <- melt(meth.geno)
        colnames(meth.geno) <- c("MK", "POS", "methGT")

        meth.error <- melt(meth.error)
        colnames(meth.error) <- c("MK", "POS", "methError")

        error.info <- merge(gab.geno, meth.geno)
        error.info <- merge(error.info, meth.error)
        error.info <- error.info[order(error.info[,2], as.character(error.info[,1])),]

        write.table(error.info, file = paste0("~{methodName}", "_updog_error_info.txt"), row.names = F, quote = F)

        if(~{methodName} == "stacks"){
           supermassa.aval <- supermassa_error(vcfR.object=vcf,
                                               onemap.object = df,
                                               vcf.par = "AD",
                                               parent1 = "P1_rg",
                                               parent2 = "P2_rg",
                                               f1="F1_rg",
                                               recovering = TRUE,
                                               mean_phred = 20,
                                               cores = 6,
                                               depths = NULL)
        } else {
          supermassa.aval <- supermassa_error(vcfR.object=vcf,
                                              onemap.object = df,
                                              vcf.par = "AD",
                                              parent1 = "P1",
                                              parent2 = "P2",
                                              f1="F1",
                                              recovering = TRUE,
                                              mean_phred = 20,
                                              cores = 6,
                                              depths = NULL)
        }
## Filters
n.mk <- supermassa.aval[[3]]
segr <- test_segregation(supermassa.aval)
distorted <- length(select_segreg(segr, distorted = T))
bins <- find_bins(supermassa.aval)
nbins <- length(bins[[1]])

filters <- data.frame("n_markers"= n.mk,
"distorted_markers"=distorted,
"redundant_markers"=n.mk-nbins)
write.table(filters, file = paste0(~{methodName},"_filters_supermassa.txt"), row.names = F, quote = F)

# Maps supermassa
twopts <- rf_2pts(supermassa.aval)
true.mks <- which(supermassa.aval[[9]] %in% tot_mks[,2])
seq.true <- make_seq(twopts, true.mks)
map.df <- map(seq.true)
map.info <- data.frame("mk.name"= colnames(supermassa.aval[[1]])[map.df[[1]]],
"pos" = supermassa.aval[[9]][map.df[[1]]],
"rf" = c(0,cumsum(haldane(map.df[[3]]))))
write.table(map.info, file = paste0(~{methodName}, "_map_supermassa.txt"), row.names = F, quote = F)


# Errors info
pos.inv <- which(supermassa.aval[[9]] %in% gab[[9]])

meth.geno <- supermassa.aval[[1]][,pos.inv]
meth.error <- supermassa.aval[[10]][,pos.inv]
meth.pos <- supermassa.aval[[9]][pos.inv]
colnames(meth.error) <- colnames(meth.geno) <- meth.pos
meth.geno <- melt(meth.geno)
colnames(meth.geno) <- c("MK", "POS", "methGT")

meth.error <- melt(meth.error)
colnames(meth.error) <- c("MK", "POS", "methError")

error.info <- merge(gab.geno, meth.geno)
error.info <- merge(error.info, meth.error)
error.info <- error.info[order(error.info[,2], as.character(error.info[,1])),]

write.table(error.info, file = paste0(~{methodName}, "_supermassa_error_info.txt"), row.names = F, quote = F)

        RSCRIPT
  >>>

  runtime{
    docker:"onemap:v1"
  } 
}

