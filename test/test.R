setwd("/home/cristiane/github/errors_workflow/cromwell-executions/F2/27bd7308-bf72-4029-9b07-d0b0cb0a2b49/")

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
  write.table(filters_tab, file = paste0(methodName,"_filters_",type.genotype,".txt"), row.names = F, quote = F)
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
  
  write.table(map.info, file = paste0(methodName, "_map_",type.genotype,".txt"), row.names = F, quote = F)
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
  
  
  if(methodName=="stacks"){
    meth.geno[,1] <- gsub("_rg","",meth.geno[,1])
    rownames(meth.error) <- gsub("_rg","",rownames(meth.error))
  }
  
  error.info <- merge(gab.geno, meth.geno)
  error.info <- error.info[order(error.info[,1], error.info[,2]),]
  error.info <- cbind(error.info, meth.error)
  
  write.table(error.info, file = paste0(methodName, "_error_info_",type.genotype,".txt"), row.names = F, quote = F)
}


tot_mks <- read.table("call-pedsim_files/execution/tot_mks.txt")
methodNames <- c("gatk")
 
for(i in 1:length(methodNames)){
  methodName <- methodNames[i]
  if(methodName == "gatk"){
    vcf <- read.vcfR("call-vcftools_filter/execution/gatk.recode.vcf")
  }
  if(methodName == "freebayes"){
    vcf <- read.vcfR("call-vcftools_filter/execution/freebayes.recode.vcf")
  }
  if(methodName == "stacks"){
    vcf <- read.vcfR("call-vcftools_filter/execution/stacks.recode.vcf")
  }
  
  if(methodName == "stacks"){
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
  
  if(methodName == "stacks"){
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
  simu <- read.vcfR("call-pedsim2vcf/execution/simu.vcf")
  gab <- onemap_read_vcfR(vcfR.object = simu,
                          cross = "f2 intercross",
                          parent1 = "P1",
                          parent2 = "P2",
                          f1 = "F1")
  
  
  
  errors_info(aval.gq, type.genotype = "GQ")
  
  # Updog
  
  if(methodName == "stacks"){
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
  if(methodName == "stacks"){
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
  
  if(methodName == "gatk"){
    polyrad.aval <- polyRAD_error(vcf="call-vcftools_filter/execution/gatk.recode.vcf", 
                                  onemap.obj = df,
                                  parent1="P1",
                                  parent2="P2",
                                  f1="F1",
                                  crosstype="f2 intercross")
  }
  if(methodName == "freebayes"){
    polyrad.aval <- polyRAD_error(vcf="call-vcftools_filter/execution/freebayes.recode.vcf", 
                                  onemap.obj = df,
                                  parent1="P1",
                                  parent2="P2",
                                  f1="F1",
                                  crosstype="f2 intercross")
  }
  if(methodName == "stacks"){
    polyrad.aval <- polyRAD_error(vcf="call-vcftools_filter/execution/stacks.recode.vcf", 
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
  
  maps(onemap.obj=polyrad.aval, type.genotype = "polyrad")
  
  ## Errors info
  errors_info(polyrad.aval, type.genotype = "polyrad")
}
