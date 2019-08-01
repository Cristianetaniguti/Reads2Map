library(supermassa4onemap)
library(onemap)
library(updog)
library(reshape2)
library(vcfR)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)

method_name <- "freebayes"
tot_mks_file <- "tot_mks.txt"
simu_vcf_file <- "simu.vcf"
vcf_file <- "freebayes.recode.vcf"

# Functions
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

create_maps_report <- function(onemap_obj, tot_mks) {
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

write_report <- function(filters_tab, out_name) {
  write.table(filters_tab, file=out_name, row.names=F, quote=F)
}

## KNOWN VARIANTS
tot_mks <- read.table(tot_mks_file)

# READING DATA FROM SIMULATED POPULATION
simu <- read.vcfR(simu_vcf_file)
gab <- onemap_read_vcfR(vcfR.object=simu,
                        cross="f2 intercross",
                        parent1="P1",
                        parent2="P2",
                        f1="F1")

## READING FINAL VCF FROM PIPELINE
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
maps_tab <- create_maps_report(df, tot_mks)
write_report(maps_tab, out_name)

# MAPS REPORT - GQ
aval.gq <- extract_depth(vcfR.object=vcf,
                         onemap.object=df,
                         vcf.par="GQ",
                         parent1="P1",
                         parent2="P2",
                         f1="F1",
                         recovering=FALSE)

aval.gq <- create_probs(df, genotypes_errors=aval.gq)

out_name <- paste0(method_name, "_map_GQ.txt")
maps_gq_tab <- create_maps_report(aval.gq, tot_mks)
write_report(maps_gq_tab, out_name)

out_name <- paste0(method_name, "_error_info_GQ.txt")
errors_tab <- create_errors_report(aval.gq, gab)
write_report(errors_tab, out_name)


# OTHER TOOLS
updog.aval <- updog_error(
  vcfR.object=vcf,
  onemap.object=df,
  vcf.par="AD",
  parent1="P1",
  parent2="P2",
  f1="F1",
  recovering=TRUE,
  mean_phred=20,
  cores=3,
  depths=NULL)

supermassa.aval <- supermassa4onemap::supermassa_error(
  vcfR.object=vcf,
  onemap.object = df,
  vcf.par = "AD",
  parent1 = "P1",
  parent2 = "P2",
  f1="F1",
  recovering = TRUE,
  mean_phred = 20,
  cores = 3,
  depths = NULL)

polyrad.aval <- polyRAD_error(
  vcf=vcf_file, 
  onemap.obj=df,
  parent1="P1",
  parent2="P2",
  f1="F1",
  crosstype="f2 intercross")

metodologies <- list(updog = updog.aval, supermassa = supermassa.aval, polyrad = polyrad.aval)
for (metodology in names(metodologies)){
  error_aval <- metodologies[[metodology]]
  ## Filters
  out_name <- paste0(method_name, "_filters_", metodology, ".txt")
  filters_tab <- create_filters_report(error_aval)
  write_report(filters_tab, out_name)
  
  ## Maps updog
  out_name <- paste0(method_name, "_map_", metodology, ".txt")
  maps_tab <- create_maps_report(error_aval, tot_mks)
  write_report(maps_tab, out_name)
  
  ## Errors info
  out_name <- paste0(method_name, "_error_", metodology, ".txt")
  errors_tab <- create_errors_report(error_aval, gab)
  write_report(errors_tab, out_name)
}

