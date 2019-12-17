# Temporary file: it will be part of empirical.wdl. I am testing separated to save time. 

library(supermassa4onemap)
library(onemap)
library(updog)
library(reshape2)
library(vcfR)
library(doParallel)
library(GUSMap)
library(ggplot2)

# Functions
source("functions.R")

  vcf_file <- "gatk_sub.vcf"
  SNPCall <- "gatk"
  cross <- "outcross"
  CountsFrom <- "vcf"
  max.cores <- 6
  parent1 <- "PT_F"
  parent2 <- "PT_M"

## READING FINAL VCF FROM PIPELINE
vcf <- read.vcfR(vcf_file)
df <- onemap_read_vcfR(vcfR.object=vcf,
                       cross= cross,
                       parent1= parent1,
                       parent2= parent2)

# removing markers with more than 75% of missing data
df <- filter_missing(df, threshold = 0.25)

# check depths
p <- create_depths_profile(onemap.obj = df, vcfR.object = vcf, parent1 = parent1, parent2 = parent2, vcf.par = "AD", 
                      recovering = FALSE, GTfrom = "vcf", alpha=0.1)
ggsave(filename = paste0(SNPCall,"_", GenoCall= "df","_", CountsFrom="vcf","_vcf_depths.png"), p)

p <- create_depths_profile(onemap.obj = df, vcfR.object = vcf, parent1 = parent1, parent2 = parent2, vcf.par = "AD", 
                           recovering = FALSE, GTfrom = "onemap", alpha=0.1)
ggsave(filename = paste0(SNPCall,"_",GenoCall="df","_",CountsFrom="vcf","_onemap_depths.png"), p)


## FILTERS REPORT
out_name <- paste0(SNPCall, "_filters_", CountsFrom,"_dfAndGQ.txt")
filters_tab <- create_filters_report(onemap_obj = df, CountsFrom, SNPCall=SNPCall, GenoCall="df")
write_report(filters_tab[[1]], out_name)


# MAPS REPORT - DF
create_map_report(input.seq = filters_tab[[2]], CountsFrom = CountsFrom, SNPCall = SNPCall, GenoCall="df")

# MAPS REPORT - GQ
aval.gq <- extract_depth(vcfR.object=vcf,
                         onemap.object=df,
                         vcf.par="GQ",
                         parent1=parent1,
                         parent2=parent2,
                         recovering=FALSE)

aval.gq <- create_probs(df, genotypes_errors=aval.gq)
filters_tab <- create_filters_report(aval.gq, CountsFrom, SNPCall=SNPCall, GenoCall="df")

create_map_report(filters_tab[[2]], CountsFrom = CountsFrom, SNPCall = SNPCall, GenoCall="GQ")

# OTHER TOOLS
## With depths from vcf

updog.aval <- updog_error(
  vcfR.object=vcf,
  onemap.object=df,
  vcf.par="AD",
  parent1=parent1,
  parent2=parent2,
  recovering=TRUE,
  mean_phred=20,
  cores=max.cores,
  depths=NULL)

supermassa.aval <- supermassa4onemap::supermassa_error(
  vcfR.object=vcf,
  onemap.object = df,
  vcf.par = "AD",
  parent1 = parent1,
  parent2 = parent2,
  recovering = TRUE,
  mean_phred = 20,
  cores = max.cores,
  depths = NULL)

polyrad.aval <- polyRAD_error(
  vcf=vcf_file,
  onemap.obj=df,
  parent1=parent1,
  parent2=parent2,
  crosstype=cross)

metodologies <- list(updog = updog.aval, supermassa = supermassa.aval, polyrad = polyrad.aval)
for (metodology in names(metodologies)){
  cat(metodology, "\n")
  error_aval <- metodologies[[metodology]]
  ## Filters
  out_name <- paste0(SNPCall, "_filters_", CountsFrom, "_",metodology, ".txt")
  filters_tab <- create_filters_report(error_aval, CountsFrom = CountsFrom, SNPCall = SNPCall, GenoCall = metodology)
  write_report(filters_tab[[1]], out_name)
  
  # check depths
  p <- create_depths_profile(onemap.obj = error_aval, vcfR.object = vcf, parent1 = parent1, parent2 = parent2, vcf.par = "AD", 
                             recovering = FALSE, GTfrom = "onemap", alpha=0.1)
  ggsave(filename = paste0(SNPCall,"_", GenoCall= metodology, "_",CountsFrom,"_onemap_depths.png"), p)
  
  ## Maps
  create_map_report(input.seq = filters_tab[[2]], CountsFrom = CountsFrom, SNPCall = SNPCall, GenoCall= metodology)
}

## Depths from bam
depths.alt <- read.table(paste0(SNPCall, "_alt_depth_bam.txt"), header = T)
depths.ref <- read.table(paste0(SNPCall, "_ref_depth_bam.txt"), header = T)

depths <- list("ref"=depths.ref, "alt"=depths.alt)

updog.aval.bam <- updog_error(
  vcfR.object=vcf,
  onemap.object=df,
  vcf.par="AD",
  parent1=parent1,
  parent2=parent2,
  recovering=TRUE,
  mean_phred=20,
  cores=max.cores,
  depths=depths)

supermassa.aval.bam <- supermassa_error(
  vcfR.object=vcf,
  onemap.object = df,
  vcf.par = "AD",
  parent1 = parent1,
  parent2 = parent2,
  recovering = TRUE,
  mean_phred = 20,
  cores = max.cores,
  depths = depths)

new.vcf <- make_vcf(vcf_file, depths, SNPCall)

polyrad.aval.bam <- polyRAD_error(
  vcf=new.vcf,
  onemap.obj=df,
  parent1=parent1,
  parent2=parent2,
  crosstype=cross)

metodologies <- list(updog = updog.aval.bam, supermassa= supermassa.aval.bam, polyrad=polyrad.aval.bam)
for (metodology in names(metodologies)){
  error_aval <- metodologies[[metodology]]
  ## Filters
  out_name <- paste0(SNPCall, "_filters_bam_", metodology, ".txt")
  filters_tab <- create_filters_report(error_aval)
  write_report(filters_tab, out_name)
  
  # check depths
  p <- create_depths_profile(onemap.obj = error_aval, vcfR.object = vcf, parent1 = parent1, parent2 = parent2, vcf.par = "AD", 
                             recovering = FALSE, GTfrom = "onemap", alpha=0.1)
  ggsave(filename = paste0(SNPCall,"_", GenoCall= metodology, "_",CountsFrom,"_onemap_depths.png"), p)
  
  ## Maps
  create_map_report(input.seq = filters_tab[[2]], CountsFrom = CountsFrom, SNPCall = SNPCall, GenoCall= metodology)
}

## Gusmap maps
out_name <- paste0(SNPCall, "_map_gusmap.txt")
map_gus <- create_gusmap_report(vcf_file, SNPCall, parent1, parent2)
write_report(map_gus, out_name)

out_name <- paste0(SNPCall, "_map_bam_gusmap.txt")
map_gus <- create_gusmap_report(new.vcf, SNPCall, parent1, parent2)
write_report(map_gus, out_name)


## Joint tables


