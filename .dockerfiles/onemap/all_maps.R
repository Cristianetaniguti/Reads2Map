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
vcf_file <- args[4]
# ~{freebayesVCF}
vcf_file <- args[5]

# print(c(method_name, tot_mks_file, simu_vcf_file, gatk_vcf_file, freebayes_vcf_file))

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

create_filters_report <- function(onemap_obj) {
    segr <- onemap::test_segregation(onemap_obj)
    distorted <- onemap::select_segreg(segr, distorted = T)
    bins <- onemap::find_bins(onemap_obj)
    total_variants <- onemap_obj[[3]]
    filters_tab <- data.frame("n_markers"= total_variants, 
                              "distorted_markers"=length(distorted), 
                              "redundant_markers"=total_variants - length(bins))
    return(filters_tab)
}

write_report <- function(filters_tab, out_name) {
    write.table(filters_tab, file=out_name, row.names=F, quote=F)
}

create_maps_report <- function(onemap_obj) {
    assign("onemap_obj", onemap_obj, envir=.GlobalEnv)
    twopts <- rf_2pts(onemap_obj)
    assign("twopts", twopts, envir=.GlobalEnv)

    true_mks <- which(onemap_obj[[9]] %in% tot_mks[,2])
    seq_true <- make_seq(twopts, true_mks)
    map_df <- map(seq_true)
    map_info <- data.frame("mk.name"= colnames(onemap_obj[[1]])[map_df[[1]]], 
                           "pos" = onemap_obj[[9]][map_df[[1]]],
                           "rf" = c(0,cumsum(haldane(map_df[[3]]))))
    return (map_info)
}

maps <- function(onemap.obj, type.genotype=NULL, method_name){
    assign("onemap.obj", onemap.obj, envir=.GlobalEnv)
    twopts <- rf_2pts(onemap.obj)
    assign("twopts", twopts, envir=.GlobalEnv)
    true.mks <- which(onemap.obj[[9]] %in% tot_mks[,2])
    seq.true <- make_seq(twopts, true.mks)
    map.df <- map(seq.true)
    map.info <- data.frame("mk.name"= colnames(onemap.obj[[1]])[map.df[[1]]], 
                           "pos" = onemap.obj[[9]][map.df[[1]]],
                           "rf" = c(0,cumsum(haldane(map.df[[3]]))))
    
    write.table(map.info, file = paste0(method_name, "_map_",type.genotype,".txt"), row.names = F, quote = F)
}

create_errors_report <- function(onemap_obj, gab) {
    pos <- which(gab[[9]] %in% onemap_obj[[9]])
    pos.inv <- which(onemap_obj[[9]] %in% gab[[9]])
    gab.pos <- gab[[9]][pos]
    gab.geno <- gab[[1]][,pos]
    colnames(gab.geno) <- gab.pos
    gab.geno <-reshape2::melt(gab.geno)
    colnames(gab.geno) <- c("MK", "POS", "gabGT")
    meth.geno <- onemap_obj[[1]][,pos.inv]
    meth.error <- onemap_obj[[11]][pos.inv + rep(c(0:(onemap_obj[[2]]-1))*onemap_obj[[3]], each=length(pos.inv)),]
    meth.pos <- onemap_obj[[9]][pos.inv]
    colnames(meth.geno) <- meth.pos
    meth.geno <- reshape2::melt(meth.geno)
    colnames(meth.geno) <- c("MK", "POS", "methGT")
    meth.error <- cbind(rownames(meth.error), meth.error)
    error.info <- merge(gab.geno, meth.geno)
    error.info <- error.info[order(error.info[,1], error.info[,2]),]
    error.info <- cbind(error.info, meth.error)
    return (error.info)
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


# READING DATA
tot_mks <- read.table(tot_mks_file)
vcf <- read.vcfR(vcf_file)
df <- onemap_read_vcfR(vcfR.object=vcf, 
                       cross="f2 intercross", 
                       parent1="P1", 
                       parent2="P2", 
                       f1="F1")

## FILTERS REPORT
out_name <- paste0(method_name, "_filters_dfAndGQ.txt")
filters_tab <- create_filters_report(df)
write_report(filters_tab, out_name)

## MAPS REPORT - DF
out_name <- paste0(method_name, "_map_df.txt")
maps_tab <- create_maps_report(df)
write_report(maps_tab, out_name)


# MAPS REPORT - GQ
aval.gq <- extract_depth(vcfR.object= vcf,
                        onemap.object = df,
                        vcf.par = "GQ",
                        parent1 = "P1",
                        parent2 = "P2",
                        f1="F1",
                        recovering = FALSE)

aval.gq <- create_probs(df, genotypes_errors= aval.gq)
out_name <- paste0(method_name, "_map_GQ.txt")
maps_gq_tab <- create_maps_report(aval.gq)
write_report(maps_gq_tab, out_name)


## READING SIMULATED VCF TO OBTAIN ERRORS
simu <- read.vcfR(simu_vcf_file)
gab <- onemap_read_vcfR(vcfR.object=simu,
                        cross="f2 intercross",
                        parent1="P1",
                        parent2="P2",
                        f1="F1")

out_name <- paste0(method_name, "_error_info_GQ.txt")
errors_tab <- create_errors_report(aval.gq, gab)
write_report(errors_tab, out_name)



# Updog
updog.aval <- updog_error(vcfR.object=vcf,
                          onemap.object = df,
                          vcf.par="AD",
                          parent1="P1",
                          parent2="P2",
                          f1="F1",
                          recovering=TRUE,
                          mean_phred=20,
                          cores=3,
                          depths=NULL)

## Filters
out_name <- paste0(method_name, "_filters_updog.txt")
filters_tab <- create_filters_report(updog.aval)
write_report(filters_tab, out_name)

## Maps updog
out_name <- paste0(method_name, "_map_updog.txt")
maps_tab <- create_maps_report(updog.aval)
write_report(maps_tab, out_name)

## Errors info
out_name <- paste0(method_name, "_error_updog.txt")
errors_tab <- create_errors_report(updog.aval, gab)
write_report(errors_tab, out_name)


# Supermassa
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

## Filters
out_name <- paste0(method_name, "_filters_supermassa.txt")
filters_tab <- create_filters_report(supermassa.aval)
write_report(filters_tab, out_name)

## Maps supermassa
out_name <- paste0(method_name, "_map_supermassa.txt")
maps_tab <- create_maps_report(supermassa.aval)
write_report(maps_tab, out_name)

## Errors info
out_name <- paste0(method_name, "_error_supermassa.txt")
errors_tab <- create_errors_report(supermassa.aval, gab)
write_report(errors_tab, out_name)

# PolyRAD
polyrad.aval <- polyRAD_error(vcf=vcf_file, 
                            onemap.obj = df,
                            parent1="P1",
                            parent2="P2",
                            f1="F1",
                            crosstype="f2 intercross")

## Filters
out_name <- paste0(method_name, "_filters_polyrad.txt")
filters_tab <- create_filters_report(polyrad.aval)
write_report(filters_tab, out_name)

## Maps supermassa
out_name <- paste0(method_name, "_map_polyrad.txt")
maps_tab <- create_maps_report(polyrad.aval)
write_report(maps_tab, out_name)

## Errors info
out_name <- paste0(method_name, "_error_polyrad.txt")
errors_tab <- create_errors_report(polyrad.aval, gab)
write_report(errors_tab, out_name)
