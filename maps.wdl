version 1.0

import "./structs/maps.wdl"

workflow maps {
  input {
    outputF2 outputf2
  }
  
  Array[String] methodNames = read_lines(outputf2.methods)  
  
  call vcftools_filter{
    input:
      gatkVCF = outputf2.gatkVCF,
      freebayesVCF = outputf2.freebayesVCF,
      stacksVCF = outputf2.stacksVCF
  }

  scatter (methodName in methodNames){
    call all_maps{
      input:
        tot_mks = outputf2.tot_mks,
        methodName = methodName,
        simu_vcf = outputf2.simu_vcf,
        gatkVCF = vcftools_filter.gatkVCF_F,
        freebayesVCF = vcftools_filter.freebayesVCF_F,
        stacksVCF = vcftools_filter.stacksVCF_F
    }
  }
}

task vcftools_filter{
  input{
    File gatkVCF
    File freebayesVCF
    File stacksVCF
  }
  output{
    File gatkVCF_F = "gatk.recode.vcf"
    File freebayesVCF_F = "freebayes.recode.vcf"
    File stacksVCF_F = "stacks.recode.vcf"
  }
  command <<<
    vcftools --vcf "~{gatkVCF}" --max-missing 0.75  --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --out gatk
    vcftools --vcf "~{freebayesVCF}" --max-missing 0.75  --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --out freebayes
    vcftools --vcf "~{stacksVCF}" --max-missing 0.75  --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --out stacks
  >>>
  runtime{
    docker: "vcftools:v1"
  }
}

task all_maps {
  input {
    File tot_mks
    String methodName
    File simu_vcf
    File gatkVCF
    File freebayesVCF
    File stacksVCF
  }
  output{
    File filters_dfAndGQ = "~{methodName}_filters_dfAndGQ.txt"
    File filters_polyrad = "~{methodName}_filters_polyrad.txt"
    File filters_supermassa = "~{methodName}_filters_supermassa.txt"
    File filters_updog = "~{methodName}_filters_updog.txt"
    File map_df = "~{methodName}_map_df.txt"
    File map_GQ = "~{methodName}_map_GQ.txt"
    File map_polyrad = "~{methodName}_map_polyrad.txt"
    File map_supermassa = "~{methodName}_map_supermassa.txt"
    File error_info_GQ = "~{methodName}_error_info_GQ.txt"
    File error_info_updog = "~{methodName}_error_info_updog.txt"
    File error_info_polyrad = "~{methodName}_error_info_polyrad.txt"
    File error_info_supermassa = "~{methodName}_error_info_supermassa.txt"
  }
  command <<<

        R --vanilla --no-save <<RSCRIPT
        # Packages
        library(supermassa4onemap)
        library(onemap)
        library(updog)
        library(reshape2)
        library(vcfR)
        library(doParallel)


        # Functions
        filters <- function(onemap.obj, type.genotype=NULL){
          n.mk <- onemap.obj[[3]]
          segr <- onemap::test_segregation(onemap.obj)
          distorted <- length(onemap::select_segreg(segr, distorted = T))
          bins <- onemap::find_bins(onemap.obj)
          nbins <- length(bins[[1]])
          
          filters_tab <- data.frame("n_markers"= n.mk, 
                                "distorted_markers"=distorted, 
                                "redundant_markers"=n.mk-nbins)
          write.table(filters_tab, file = paste0("~{methodName}","_filters_",type.genotype,".txt"), row.names = F, quote = F)
        }

        maps <- function(onemap.obj, type.genotype=NULL){
          assign("onemap.obj", onemap.obj, envir = .GlobalEnv)
          twopts <- rf_2pts(onemap.obj)
          assign("twopts", twopts, envir = .GlobalEnv)
          true.mks <- which(onemap.obj[[9]] %in% tot_mks[,2])
          seq.true <- make_seq(twopts, true.mks)
          map.df <- map(seq.true)
          map.info <- data.frame("mk.name"= colnames(onemap.obj[[1]])[map.df[[1]]], 
                                "pos" = onemap.obj[[9]][map.df[[1]]],
                                "rf" = c(0,cumsum(haldane(map.df[[3]]))))
          
          write.table(map.info, file = paste0("~{methodName}", "_map_",type.genotype,".txt"), row.names = F, quote = F)
        }

        errors_info <- function(onemap.obj=NULL, type.genotype=NULL){
          pos <- which(gab[[9]] %in% onemap.obj[[9]])
          pos.inv <- which(onemap.obj[[9]] %in% gab[[9]])
          
          gab.pos <- gab[[9]][pos]
          gab.geno <- gab[[1]][,pos]
          colnames(gab.geno) <- gab.pos
          gab.geno <-reshape2::melt(gab.geno)
          colnames(gab.geno) <- c("MK", "POS", "gabGT")
          
          meth.geno <- onemap.obj[[1]][,pos.inv]
          meth.error <- onemap.obj[[11]][pos.inv + rep(c(0:(onemap.obj[[2]]-1))*onemap.obj[[3]], each=length(pos.inv)),]
          meth.pos <- onemap.obj[[9]][pos.inv]
          colnames(meth.geno) <- meth.pos
          meth.geno <- reshape2::melt(meth.geno)
          colnames(meth.geno) <- c("MK", "POS", "methGT")
          meth.error <- cbind(rownames(meth.error), meth.error)
          
          
          if("~{methodName}"=="stacks"){
            meth.geno[,1] <- gsub("_rg","",meth.geno[,1])
            rownames(meth.error) <- gsub("_rg","",rownames(meth.error))
          }
          
          error.info <- merge(gab.geno, meth.geno)
          error.info <- error.info[order(error.info[,1], error.info[,2]),]
          error.info <- cbind(error.info, meth.error)
          
          write.table(error.info, file = paste0("~{methodName}", "_error_info_",type.genotype,".txt"), row.names = F, quote = F)
        }


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

        filters(onemap.obj = df, type.genotype = "dfAndGQ")

        ## Maps default

        maps(onemap.obj = df, type.genotype = "df")

        # GQ

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

        aval.gq <- create_probs(df, genotypes_errors= aval.gq)

        ## Maps GQ

        maps(aval.gq, type.genotype = "GQ")

        ## Errors info tab
        simu <- read.vcfR("~{simu_vcf}")
        gab <- onemap_read_vcfR(vcfR.object = simu,
                                cross = "f2 intercross",
                                parent1 = "P1",
                                parent2 = "P2",
                                f1 = "F1")



        errors_info(aval.gq, type.genotype = "GQ")

        # Updog

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

        filters(updog.aval, type.genotype = "updog")

        ## Maps updog

        maps(updog.aval, type.genotype = "updog")

        ## Errors info
        errors_info(updog.aval, type.genotype = "updog")

        # Supermassa
        if("~{methodName}" == "stacks"){
          supermassa.aval <- supermassa4onemap::supermassa_error(vcfR.object=vcf,
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
          supermassa.aval <- supermassa4onemap::supermassa_error(vcfR.object=vcf,
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
        filters(supermassa.aval, type.genotype="supermassa")

        ## Maps supermassa

        maps(supermassa.aval, type.genotype = "supermassa")

        ## Errors info
        errors_info(supermassa.aval, type.genotype = "supermassa")

        # PolyRAD

        if("~{methodName}" == "gatk"){
          polyrad.aval <- polyRAD_error(vcf="~{gatkVCF}", 
                                        onemap.obj = df,
                                        parent1="P1",
                                        parent2="P2",
                                        f1="F1",
                                        crosstype="f2 intercross")
        }
        if("~{methodName}" == "freebayes"){
          polyrad.aval <- polyRAD_error(vcf="~{freebayesVCF}", 
                                        onemap.obj = df,
                                        parent1="P1",
                                        parent2="P2",
                                        f1="F1",
                                        crosstype="f2 intercross")
        }
        if("~{methodName}" == "stacks"){
          polyrad.aval <- polyRAD_error(vcf="~{stacksVCF}", 
                                        onemap.obj = df,
                                        parent1="P1_rg",
                                        parent2="P2_rg",
                                        f1="F1_rg",
                                        crosstype="f2 intercross",
                                        tech.issue=F)
        }

        ## Filters
        filters(polyrad.aval, type.genotype="polyrad")

        ## Maps supermassa

        maps(polyrad.aval, type.genotype = "polyrad")

        ## Errors info
        errors_info(polyrad.aval, type.genotype = "polyrad")

        RSCRIPT
  >>>

  runtime{
    docker:"onemap:v1"
  } 
}

