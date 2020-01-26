library(R.utils)

tsv_files <- read.table("counts", stringsAsFactors = F)
tsv_files <- tsv_files$V1
freebayes_counts <- tsv_files[grep("freebayes", tsv_files)]
gatk_counts <- tsv_files[grep("gatk", tsv_files)]
counts <- list(freebayes_counts, gatk_counts)
names <- sapply(strsplit(sapply(strsplit(gatk_counts, "/"), "[",10), "_gatk"), "[", 1)
  
  
methods <- c("freebayes", "gatk")

# BamCounts4OneMap
for(method in 1:2){
  
  file.counts <- read.table(counts[[method]][1], skip = 1448, header=T, stringsAsFactors = F)
  
  ref_depth_matrix2 <- alt_depth_matrix2  <- matrix(NA, nrow = dim(file.counts)[1], ncol = length(counts[[method]]))
  
  for(j in 1:length(counts[[method]])){
    ## From picard tool
    
    file.counts <- read.table(counts[[method]][j], skip = 1448, header=T, stringsAsFactors = F)
    
    ref_depth_matrix2[,j] <- file.counts[,3]
    alt_depth_matrix2[,j] <- file.counts[,4]
    
    if(j == 1){
      ref_allele <- file.counts[,5]
      alt_allele <- file.counts[,6]
    } else {
      idx.ref <- which(ref_allele == "N")
      idx.alt <- which(alt_allele == "N")
      if(length(idx.ref)!=0){
        ref_allele[idx.ref] <- file.counts[idx.ref,5]
      }
      if(length(idx.alt)!=0){
        alt_allele[idx.alt] <- file.counts[idx.alt,6]
      }
    }
    
  }
  
  rownames(ref_depth_matrix2) <- rownames(alt_depth_matrix2) <- paste0(file.counts[,1],"_", file.counts[,2])
  colnames(ref_depth_matrix2) <- colnames(alt_depth_matrix2) <- names
  
  alleles <- data.frame(file.counts[,1],file.counts[,2], ref_allele, alt_allele)
  write.table(alleles, file = paste0(methods[method],"_example4ref_alt_alleles.txt"), col.names = F, row.names = F)
  
  write.table(ref_depth_matrix2, file = paste0(methods[method],"_ref_depth_bam.txt"), quote=F, row.names=T, sep="\t", col.names=T)
  write.table(alt_depth_matrix2, file = paste0(methods[method],"_alt_depth_bam.txt"), quote=F, row.names=T, sep="\t", col.names=T)
}


##########################################################################
library(supermassa4onemap)
library(onemap)
library(updog)
library(reshape2)
library(vcfR)
library(doParallel)
library(GUSMap)

method_name <- c("freebayes", "gatk")
vcf_files <- c("freebayes.recode.vcf", "gatk.recode.vcf")
vcf_file <- "gatk.recode.vcf"
cross <- "outcross"

# Functions
parmap <- function(input.seq=NULL, cores=3, overlap=4, tol=10E-5, avoid_link_errors = TRUE){
  twopts <- input.seq$twopt
  
  interv <- length(input.seq$seq.num)/cores
  
  seqs <- 1:length(input.seq$seq.num)
  
  list_idx <- list(1:cores)
  
  init <- 1
  end <- interv
  for(i in 1:cores){
    if(i != cores){
      list_idx[[i]] <- init:(end+(overlap-1))
      init <- end
      end <- end + interv
    } else{
      list_idx[[i]] <- init:end
    }
  }
  
  list_seq <- lapply(list_idx, function(x) make_seq(twopts, input.seq$seq.num[x]))
  

  clust <- makeCluster(cores)
  clusterExport(clust, c("avoid_unlinked", "twopts"))
  if(avoid_link_errors){
    new.maps <- parLapply(clust, list_seq, function(x) avoid_unlinked(x, tol=tol))      
  } else {
    new.maps <- parLapply(clust, list_seq, function(x) onemap::map(x, tol=tol))
  }

  stopCluster(clust)
  
  joint.map <- new.maps[[1]]
  new.seq.num <- new.seq.rf <- new.seq.phases <- vector()
  diff1 <- vector()
  for(i in 1:(length(new.maps)-1)){
    idx.end <- (length(new.maps[[i]]$seq.num) - overlap+1):(length(new.maps[[i]]$seq.num))
    new.seq.num <- c(new.seq.num, new.maps[[i]]$seq.num[-idx.end])
    new.seq.rf <- c(new.seq.rf, new.maps[[i]]$seq.rf[-idx.end[-length(idx.end)]])
    new.seq.phases <- c(new.seq.phases, new.maps[[i]]$seq.phases[-idx.end[-length(idx.end)]])
    
    if(i == (length(new.maps)-1)){
      new.seq.num <- c(new.seq.num, new.maps[[i+1]]$seq.num)
      new.seq.rf <- c(new.seq.rf, new.maps[[i+1]]$seq.rf)
      new.seq.phases <- c(new.seq.phases, new.maps[[i+1]]$seq.phases)
    }
    
    end <- new.maps[[i]]$seq.rf[idx.end[-length(idx.end)]]
    init <- new.maps[[i+1]]$seq.rf[1:overlap-1]
    diff1 <- c(diff1,end - init)
  }
  
  cat("The overlap markers have mean ", mean(diff1), " of  recombination fraction diff1erences, and variance of ", var(diff1), "\n")
  
  joint.map$seq.num <- new.seq.num
  joint.map$seq.phases <- new.seq.phases
  joint.map$seq.rf <- new.seq.rf
  
  return(list(diff1,joint.map))
}

avoid_unlinked <- function(input.seq, tol){
    map_df <- onemap::map(input.seq, mds.seq = T)
    while(class(map_df) == "integer"){
        seq_true <- onemap::make_seq(twopts, map_df)
        map_df <- onemap::map(input.seq = seq_true, mds.seq = T)
    }
    return(map_df)
}


create_filters_report <- function(onemap_obj) {
  bins <- onemap::find_bins(onemap_obj)
  onemap_bins <- create_data_bins(onemap_obj, bins)
  segr <- onemap::test_segregation(onemap_bins)
  distorted <- onemap::select_segreg(segr, distorted = T)
  no_distorted <- onemap::select_segreg(segr, distorted = F, numbers = T)
  twopts <- rf_2pts(onemap_bins)
  seq1 <- make_seq(twopts, no_distorted)
  total_variants <- onemap_obj[[3]]
  filters_tab <- data.frame("n_markers"= total_variants,
                            "distorted_markers"=length(distorted),
                            "redundant_markers"=total_variants - length(bins[[1]]))
  return(list(filters_tab, seq1))
}

write_report <- function(filters_tab, out_name) {
  write.table(filters_tab, file=out_name, row.names=F, quote=F)
}


make_vcf <- function(vcf.old, depths, method){
  # The input od polyRAD need to be a VCF, then this part takes the allele depth from "depths" and put at AD field of input vcf
  idx <- system(paste0("grep -in 'CHROM' ", vcf.old), intern = T) # This part only works in linux OS
  idx.i <- strsplit(idx, split = ":")[[1]][1]
  seed <- sample(1:10000, 1)
  system(paste0("head -n ", idx.i," ", vcf.old, " > head.",seed))
  
  vcf.tab <- read.table(vcf.old, stringsAsFactors = F)
  
  if(all(rownames(depths[[1]]) == paste0(vcf.tab[,1], "_", vcf.tab[,2]))){
    
    vcf.init <- vcf.tab[,1:8]
    AD.colum <- rep("AD", dim(vcf.init)[1])
    vcf.init <- cbind(vcf.init, AD.colum)
    
    rs <- rownames(depths[[1]])
    vcf.init[,3] <- rs
  } else {
    temp.tab <- read.table(paste0(method,"_example4ref_alt_alleles.txt"))
    vcf.init <- cbind(temp.tab[,1:2],paste0(temp.tab[,1], "_", temp.tab[,2]), temp.tab[,3:4],
                      rep(".", dim(temp.tab)[1]),rep(".", dim(temp.tab)[1]), rep(".", dim(temp.tab)[1]), rep("AD",dim(temp.tab)[1]))
  }
  
  ind.n <- colnames(depths[[1]]) # The names came in different order
  
  header <- strsplit(idx, split = "\t")[[1]]
  ind.vcf <- header[10:length(header)]
  ind.n <- factor(ind.n, levels = ind.vcf)
  
  depths[[1]] <- depths[[1]][,order(ind.n)]
  depths[[2]] <- depths[[2]][,order(ind.n)]
  
  comb.depth <- matrix(paste0(as.matrix(depths[[1]]), ",", as.matrix(depths[[2]])), ncol = ncol(depths[[2]]))
  colnames(comb.depth) <- ind.vcf
  #hmc.file <- cbind(rs, comb.depth)
  
  vcf.body <- cbind(vcf.init, comb.depth)
  
  write.table(vcf.body, file = paste0("temp.body.", seed), quote = FALSE, sep = "\t", row.names = FALSE, col.names = F)
  
  system(paste0("cat head.",seed," temp.body.",seed," > temp.",seed,".vcf"))
  return(paste0("temp.",seed, ".vcf"))
}

create_gusmap_report <- function(vcf_file){
  ## Maps with gusmap
  RAfile <- VCFtoRA(vcf_file, makePed = T)
  
  ped.file <- read.csv(paste0(method_name,"_ped.csv"))
  ID <- c(1:(dim(ped.file)[1]))
  idx.P1 <- which(as.character(ped.file$SampleID) == "PT_F")
  idx.P2 <- which(as.character(ped.file$SampleID) == "PT_M")
  
  mother <- c(rep(idx.P1,dim(ped.file)[1]))
  father <- c(rep(idx.P2,dim(ped.file)[1]))
  fam <- c(rep("F1", dim(ped.file)[1]))
  idx.P1 <- which(as.character(ped.file$SampleID) == "PT_F")
  idx.P2 <- which(as.character(ped.file$SampleID) == "PT_M")
  mother[c(idx.P1, idx.P2)] <- ""
  father[c(idx.P1, idx.P2)] <- ""
  fam[c(idx.P1, idx.P2)] <- c("", "")
  ped.file$IndividualID <- ID
  ped.file$Mother <- father # Inverted
  ped.file$Father <- mother
  ped.file$Family <- fam
  
  write.csv(ped.file, file = "ped.file.csv")
  
  RAdata <- readRA(paste0(method_name,".recode.vcf.ra.tab"), pedfile = "ped.file.csv", 
                   filter = list(MAF=0.05, MISS=0.75, BIN=0, DEPTH=0, PVALUE=0.01), sampthres = 0)
  
  mydata <- makeFS(RAobj = RAdata, pedfile = "ped.file.csv", 
                   filter = list(MAF = 0.05, MISS = 0.75,
                                 BIN = 0, DEPTH = 0, PVALUE = 0.01))

  
  pos <- mydata$.__enclos_env__$private$pos
  depth_Ref_m <- mydata$.__enclos_env__$private$ref[[1]]
  depth_Alt_m <- mydata$.__enclos_env__$private$alt[[1]]
  
  depth_Ref <- list(depth_Ref_m)
  depth_Alt <- list(depth_Alt_m)
  config <- mydata$.__enclos_env__$private$config[[1]]
  
  phases.gus <- GUSMap:::infer_OPGP_FS(depth_Ref_m, depth_Alt_m, 
                                       config, epsilon=0.001, reltol=1e-3)
  
  rf_est <- GUSMap:::rf_est_FS(init_r = 0.01, ep = 0.001, 
                               ref = depth_Ref, 
                               alt = depth_Alt, 
                               OPGP=list(as.integer(phases.gus)),
                               nThreads = 1)
  
  rf_est$rf[which(rf_est$rf > 0.5)] <- 0.4999999
  dist.gus <- c(0,cumsum(rf_est$rf))
  phases.gus[which(phases.gus == 1 | phases.gus == 4)] <- 17
  phases.gus[which(phases.gus == 2 | phases.gus == 3)] <- 18
  phases.gus[which(phases.gus == 5 | phases.gus == 8)] <- 19
  phases.gus[which(phases.gus == 6 | phases.gus == 7)] <- 20
  phases.gus[which(phases.gus == 9 | phases.gus == 12)] <- 21
  phases.gus[which(phases.gus == 10 | phases.gus == 11)] <- 22
  phases.gus[which(phases.gus == 13 | phases.gus == 16)] <- 23
  phases.gus[which(phases.gus == 14 | phases.gus == 15)] <- 24

  config[which(config==1)] <- "B3.7"
  config[which(config==2 | config==3)] <- "D1.10"
  config[which(config==4 | config==5)] <- "D2.15"
  
  map_info <- data.frame("mk.name"= mydata$.__enclos_env__$private$SNP_Names,
                         "pos" = pos,
                         "rf" = dist.gus,
                         "type"= config,
                         "est.phases"= phases.gus)
  
  return(map_info)
}

vcf_file <- "gatk.recode.vcf"
method_name <- "gatk"
## READING FINAL VCF FROM PIPELINE
vcf <- read.vcfR(vcf_file)
df <- onemap_read_vcfR(vcfR.object=vcf,
                       cross= cross,
                       parent1="PT_F",
                       parent2="PT_M")


filter_missing <- function(df=NULL, threshold= 0.25){
  perc.mis <- apply(df$geno, 2, function(x) sum(x == 0)/length(x))
  idx <- which(!perc.mis > threshold)
  
  new.df <- df
  new.df$geno <- df$geno[,idx]
  new.df$n.mar <- length(idx)
  new.df$segr.type <- df$segr.type[idx]
  new.df$segr.type.num <- df$segr.type.num[idx]
  new.df$CHROM <- df$CHROM[idx]
  new.df$POS <- df$POS[idx]
  new.df$error <- df$error[idx + rep(c(0:(df$n.ind-1))*df$n.mar, each=length(idx)),]
  cat("Number of markers removed from the onemap object: ", length(which(perc.mis > threshold)))
  return(new.df)
}

# removing markers with more than 75% of missing data
df <- filter_missing(df, threshold = 0.25)

## FILTERS REPORT
out_name <- paste0(method_name, "_filters_dfAndGQ.txt")
filters_tab <- create_filters_report(df)
write_report(filters_tab[[1]], out_name)

test_seq <- make_seq(filters_tab[[2]]$twopt, filters_tab[[2]]$seq.num[1:150])

# Parallel from BatchMap
batch_size <- pick_batch_sizes(filters_tab[[2]],
                               size = 60,
                               overlap = 25,
                               around = 10)

time_batchmap <- system.time(map_out <- map_overlapping_batches(input.seq = filters_tab[[2]], 
                                                           size = batch_size, 
                                                           phase_cores = 4, 
                                                           overlap = 25))



# parmap
max.cores <- detectCores() -4 # Change here to the maximun cores available for the process
n.mk <- length(filters_tab[[2]]$seq.num)
min.group.size <- 60

while(n.mk/max.cores < min.group.size){
  max.cores <- max.cores - 1
}

out_name <- paste0(method_name, "_map_df.RData")
time_parmap <- system.time(map_out <- parmap(input.seq = filters_tab[[2]], cores = max.cores, overlap = 10))
save(map_out, out_name)

time_tot <- system.time(
  map_tot <- map(filters_tab[[2]], rm_unlinked = T)
while(class(map_tot) == "interger"){
  map_tot <- map(filters_tab[[2]], rm_unlinked = T)
})


# MAPS REPORT - GQ
aval.gq <- extract_depth(vcfR.object=vcf,
                         onemap.object=df,
                         vcf.par="GQ",
                         parent1="PT_F",
                         parent2="PT_M",
                         f1 = f1,
                         recovering=FALSE)

aval.gq <- create_probs(df, genotypes_errors=aval.gq)
filters_tab <- create_filters_report(aval.gq)

out_name <- paste0(method_name, "_map_GQ.RData")
time_parmap <- system.time(map_out <- parmap(input.seq = filters_tab[[2]], cores = max.cores, overlap = 10))
save(map_out, out_name)

# OTHER TOOLS
## With depths from vcf

updog.aval <- updog_error(
  vcfR.object=vcf,
  onemap.object=df,
  vcf.par="AD",
  parent1="PT_F",
  parent2="PT_M",
  f1 = f1,
  recovering=TRUE,
  mean_phred=20,
  cores=max.cores,
  depths=NULL)

supermassa.aval <- supermassa4onemap::supermassa_error(
  vcfR.object=vcf,
  onemap.object = df,
  vcf.par = "AD",
  parent1 = "PT_F",
  parent2 = "PT_M",
  f1 = f1,
  recovering = TRUE,
  mean_phred = 20,
  cores = max.cores,
  depths = NULL)

polyrad.aval <- polyRAD_error(
  vcf=vcf_file,
  onemap.obj=df,
  parent1="PT_F",
  parent2="PT_M",
  f1 = f1,
  crosstype=cross)

metodologies <- list(updog = updog.aval, supermassa = supermassa.aval, polyrad = polyrad.aval)
for (metodology in names(metodologies)){
  error_aval <- metodologies[[metodology]]
  ## Filters
  out_name <- paste0(method_name, "_filters_", metodology, ".txt")
  filters_tab <- create_filters_report(error_aval)
  write_report(filters_tab[[1]], out_name)
  
  ## Maps
  out_name <- paste0(method_name, "_map_", metodology, ".RData")
  map_out <- parmap(input.seq = filters_tab[[2]], cores = max.cores, overlap = 10)
  save(map_out, out_name)
}

## Depths from bam
depths.alt <- read.table(paste0(method_name, "_alt_depth_bam.txt"), header = T)
depths.ref <- read.table(paste0(method_name, "_ref_depth_bam.txt"), header = T)

depths <- list("ref"=depths.ref, "alt"=depths.alt)

updog.aval.bam <- updog_error(
  vcfR.object=vcf,
  onemap.object=df,
  vcf.par="AD",
  parent1="PT_F",
  parent2="PT_M",
  f1 = f1,
  recovering=TRUE,
  mean_phred=20,
  cores=max.cores,
  depths=depths)

supermassa.aval.bam <- supermassa_error(
  vcfR.object=vcf,
  onemap.object = df,
  vcf.par = "AD",
  parent1 = "PT_F",
  parent2 = "PT_M",
  f1 = f1,
  recovering = TRUE,
  mean_phred = 20,
  cores = max.cores,
  depths = depths)

new.vcf <- make_vcf(vcf_file, depths, method_name)

polyrad.aval.bam <- polyRAD_error(
  vcf=new.vcf,
  onemap.obj=df,
  parent1="PT_F",
  parent2="PT_M",
  f1 = f1,
  crosstype=cross)

metodologies <- list(updog = updog.aval.bam, supermassa= supermassa.aval.bam, polyrad=polyrad.aval.bam)
for (metodology in names(metodologies)){
  error_aval <- metodologies[[metodology]]
  ## Filters
  out_name <- paste0(method_name, "_filters_bam_", metodology, ".txt")
  filters_tab <- create_filters_report(error_aval)
  write_report(filters_tab, out_name)

  ## Maps
  out_name <- paste0(method_name, "_map_bam_", metodology, ".RData")
  map_out <- parmap(input.seq = filters_tab[[2]], cores = max.cores, overlap = 10)
  save(map_out, out_name)
}

## Gusmap maps
out_name <- paste0(method_name, "_map_gusmap.txt")
map_gus <- create_gusmap_report(vcf_file)
write_report(map_gus, out_name)

out_name <- paste0(method_name, "_map_bam_gusmap.txt")
map_gus <- create_gusmap_report(new.vcf)
write_report(map_gus, out_name)
