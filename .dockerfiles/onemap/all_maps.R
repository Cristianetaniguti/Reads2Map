library(supermassa4onemap)
library(onemap)
library(updog)
library(reshape2)
library(vcfR)
library(doParallel)



args = commandArgs(trailingOnly=TRUE)

# method_name <- "~{methodName}"
method_name <- args[1]
# "~{tot_mks}"
tot_mks_file <- args[2]
# ~{simu_vcf}
simu_vcf_file <- args[3]
# ~{gatkVCF}
gatk_vcf_file <- args[4]
# ~{freebayesVCF}
freebayes_vcf_file <- args[5]

print(c(method_name, tot_mks_file, simu_vcf_file, gatk_vcf_file, freebayes_vcf_file))

# Functions
filters <- function(onemap.obj, type.genotype=NULL, method_name){
    n.mk <- onemap.obj[[3]]
    segr <- onemap::test_segregation(onemap.obj)
    distorted <- length(onemap::select_segreg(segr, distorted = T))
    bins <- onemap::find_bins(onemap.obj)
    nbins <- length(bins[[1]])
    
    filters_tab <- data.frame("n_markers"= n.mk, 
                        "distorted_markers"=distorted, 
                        "redundant_markers"=n.mk-nbins)
    write.table(filters_tab, file = paste0(method_name, "_filters_", type.genotype,".txt"), row.names = F, quote = F)
}

maps <- function(onemap.obj, type.genotype=NULL, method_name){
    assign("onemap.obj", onemap.obj, envir = .GlobalEnv)
    twopts <- rf_2pts(onemap.obj)
    assign("twopts", twopts, envir = .GlobalEnv)
    true.mks <- which(onemap.obj[[9]] %in% tot_mks[,2])
    seq.true <- make_seq(twopts, true.mks)
    map.df <- map(seq.true)
    map.info <- data.frame("mk.name"= colnames(onemap.obj[[1]])[map.df[[1]]], 
                        "pos" = onemap.obj[[9]][map.df[[1]]],
                        "rf" = c(0,cumsum(haldane(map.df[[3]]))))
    
    write.table(map.info, file = paste0(method_name, "_map_",type.genotype,".txt"), row.names = F, quote = F)
}

errors_info <- function(onemap.obj=NULL, type.genotype=NULL, method_name=NULL){
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
    
    
    if(method_name=="stacks"){
        meth.geno[,1] <- gsub("_rg","",meth.geno[,1])
        rownames(meth.error) <- gsub("_rg","",rownames(meth.error))
    }
    
    error.info <- merge(gab.geno, meth.geno)
    error.info <- error.info[order(error.info[,1], error.info[,2]),]
    error.info <- cbind(error.info, meth.error)
    
    write.table(error.info, file = paste0(method_name, "_error_info_",type.genotype,".txt"), row.names = F, quote = F)
}


errors_info2 <- function(onemap.obj=NULL, type.genotype=NULL, method_name=NULL){
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
    
    
    if(method_name=="stacks"){
        meth.geno[,1] <- gsub("_rg","",meth.geno[,1])
        rownames(meth.error) <- gsub("_rg","",rownames(meth.error))
    }
    
    error.info <- merge(gab.geno, meth.geno)
    error.info <- error.info[order(error.info[,1], error.info[,2]),]
    error.info <- cbind(error.info, meth.error)
    
    write.table(error.info, file = paste0(method_name, "_error_info_",type.genotype,".txt"), row.names = F, quote = F)
}

tot_mks <- read.table(tot_mks_file)

if(method_name == "gatk"){
    vcf <- read.vcfR(gatk_vcf_file)
}
if(method_name == "freebayes"){
    vcf <- read.vcfR(freebayes_vcf_file)
}

if(method_name == "stacks"){
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
filters(onemap.obj=df, type.genotype="dfAndGQ", method_name=method_name)

## Maps default
maps(onemap.obj = df, type.genotype = "df", method_name=method_name)

# GQ
if(method_name == "stacks"){
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
maps(aval.gq, type.genotype = "GQ", method_name=method_name)

## Errors info tab
simu <- read.vcfR(simu_vcf_file)
gab <- onemap_read_vcfR(vcfR.object = simu,
                        cross = "f2 intercross",
                        parent1 = "P1",
                        parent2 = "P2",
                        f1 = "F1")



errors_info(aval.gq, type.genotype = "GQ", method_name=method_name)

# Updog
if(method_name == "stacks"){
    updog.aval <- updog_error(vcfR.object=vcf,
                            onemap.object = df,
                            vcf.par = "AD",
                            parent1 = "P1_rg",
                            parent2 = "P2_rg",
                            f1="F1_rg",
                            recovering = TRUE,
                            mean_phred = 20,
                            cores = 3,
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
                            cores = 3,
                            depths = NULL)
}

## Filters        
filters(updog.aval, type.genotype = "updog", method_name=method_name)

## Maps updog
maps(updog.aval, type.genotype = "updog", method_name=method_name)

## Errors info
errors_info(updog.aval, type.genotype = "updog", method_name=method_name)

# Supermassa
if(method_name == "stacks"){
    supermassa.aval <- supermassa4onemap::supermassa_error(vcfR.object=vcf,
                                        onemap.object = df,
                                        vcf.par = "AD",
                                        parent1 = "P1_rg",
                                        parent2 = "P2_rg",
                                        f1="F1_rg",
                                        recovering = TRUE,
                                        mean_phred = 20,
                                        cores = 3,
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
                                        cores = 3,
                                        depths = NULL)
}

## Filters
filters(supermassa.aval, type.genotype="supermassa", method_name=method_name)

## Maps supermassa

maps(supermassa.aval, type.genotype = "supermassa", method_name=method_name)

## Errors info
errors_info(supermassa.aval, type.genotype = "supermassa", method_name=method_name)

# PolyRAD

if(method_name == "gatk"){
    polyrad.aval <- polyRAD_error(vcf=gatk_vcf_file, 
                                onemap.obj = df,
                                parent1="P1",
                                parent2="P2",
                                f1="F1",
                                crosstype="f2 intercross")
}
if(method_name == "freebayes"){
    polyrad.aval <- polyRAD_error(vcf=freebayes_vcf_file, 
                                onemap.obj = df,
                                parent1="P1",
                                parent2="P2",
                                f1="F1",
                                crosstype="f2 intercross")
}

## Filters
filters(polyrad.aval, type.genotype="polyrad", method_name=method_name)

## Maps supermassa
maps(polyrad.aval, type.genotype = "polyrad", method_name=method_name)

## Errors info
errors_info(polyrad.aval, type.genotype = "polyrad", method_name=method_name)
