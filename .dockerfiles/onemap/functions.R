# Functions
parmap <- function(input.seq=NULL, cores=3, overlap=4, tol=10E-5){
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
  
  no_cores <- detectCores()
  
  clust <- makeCluster(no_cores)
  
  new.maps <- parLapply(clust, list_seq, function(x) onemap::map(x, tol=tol))
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

phaseToOPGP_OM <- function(x){
  ## code from here taken from the onemap function print.sequence()
  link.phases <- matrix(NA, length(x$seq.num), 2)
  link.phases[1, ] <- rep(1, 2)
  for (i in 1:length(x$seq.phases)) {
    switch(EXPR = x$seq.phases[i],
           link.phases[i + 1, ] <- link.phases[i, ] * c(1, 1),
           link.phases[i +  1, ] <- link.phases[i, ] * c(1, -1),
           link.phases[i + 1, ] <- link.phases[i, ] * c(-1, 1),
           link.phases[i + 1, ] <- link.phases[i, ] * c(-1, -1))
  }
  if (class(x$data.name)[2] == "outcross") {
    link.phases <- apply(link.phases, 1, function(x) paste(as.character(x), collapse = "."))
    parents <- matrix("", length(x$seq.num), 4)
    for (i in 1:length(x$seq.num)) 
      parents[i, ] <- onemap:::return_geno(x$data.name$segr.type[x$seq.num[i]], link.phases[i])
    ## Our code below
    #transpose the parents and set to baseline
    parents[which(parents == 'a')] <-'A'
    parents[which(parents == 'b')] <- 'B'
    
    parents = t(parents)
    if(parents[1,which(apply(parents[1:2,],2,function(x) !(all(x=='A'))))[1]] == 'B')
      parents[1:2,] <- parents[2:1,]
    if(parents[3,which(apply(parents[3:4,],2,function(x) !(all(x=='A'))))[1]] == 'B')
      parents[3:4,] <- parents[4:3,]
    
    phases <- GUSMap:::parHapToOPGP(parents)
    phases[which(phases == 1 | phases == 4)] <- 17 
    phases[which(phases == 2 | phases == 3)] <- 18
    phases[which(phases == 5 | phases == 8)] <- 19
    phases[which(phases == 6 | phases == 7)] <- 20
    phases[which(phases == 9 | phases == 12)] <- 21
    phases[which(phases == 10 | phases == 11)] <- 22
    phases[which(phases == 13 | phases == 16)] <- 23 
    phases[which(phases == 14 | phases == 15)] <- 24 
    
    ## Now from the parental haplotypes, determine the OPGPs
    return(phases)
  }
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

create_maps_report <- function(input.seq, tot_mks) {
  
  true_mks <- input.seq$seq.num[which(input.seq$data.name$POS[input.seq$seq.num] %in% tot_mks[,2])]
  seq_true <- make_seq(input.seq$twopt, true_mks) # only true markers are mapped
  
  max.cores <- 10 #detectCores() # Change here to the maximun cores available for the process
  n.mk <- length(input.seq$seq.num)
  min.group.size <- 60
  
  while(n.mk/max.cores < min.group.size){
    max.cores <- max.cores - 1
  }
  if(max.cores == 0) {
    map_df <- map(seq_true)
  } else {
    map_df <- parmap(input.seq = seq_true, cores = max.cores, overlap = 10)
    map_df <- map_df[[2]]
  }
  
  phases <- phaseToOPGP_OM(map_df)
  real_phase <- real_phases[which(real_phases[,1] %in% input.seq$data.name$POS[map_df[[1]]]),2]
  map_info <- data.frame("mk.name"= colnames(input.seq$data.name$geno)[map_df[[1]]],
                         "pos" = input.seq$data.name$POS[map_df[[1]]],
                         "rf" = c(0,cumsum(haldane(map_df[[3]]))),
                         "type"= input.seq$data.name$segr.type[map_df[[1]]],
                         "est.phases"= phases,
                         "real.phases"= real_phase)
  return (map_info)
}


create_gusmap_report <- function(vcf_file){
  ## Maps with gusmap
  RAfile <- VCFtoRA(vcf_file, makePed = T)
  
  ped.file <- read.csv(paste0(method_name,"_ped.csv"))
  ID <- c(1:(dim(ped.file)[1]))
  idx.P1 <- which(as.character(ped.file$SampleID) == "P1")
  idx.P2 <- which(as.character(ped.file$SampleID) == "P2")
  
  mother <- c(rep(idx.P1,dim(ped.file)[1]))
  father <- c(rep(idx.P2,dim(ped.file)[1]))
  fam <- c(rep("F1", dim(ped.file)[1]))
  idx.P1 <- which(as.character(ped.file$SampleID) == "P1")
  idx.P2 <- which(as.character(ped.file$SampleID) == "P2")
  mother[c(idx.P1, idx.P2)] <- ""
  father[c(idx.P1, idx.P2)] <- ""
  fam[c(idx.P1, idx.P2)] <- c("", "")
  ped.file$IndividualID <- ID
  ped.file$Mother <- father # Inverted
  ped.file$Father <- mother
  ped.file$Family <- fam
  
  write.csv(ped.file, file = "ped.file.csv")
  
  RAdata <- readRA(paste0(method_name,".recode.vcf.ra.tab"), pedfile = "ped.file.csv", 
                   filter = list(MAF=0.05, MISS=1, BIN=0, DEPTH=0, PVALUE=0), sampthres = 0)
  
  mydata <- makeFS(RAobj = RAdata, pedfile = "ped.file.csv", 
                   filter = list(MAF = 0.05, MISS = 1,
                                 BIN = 0, DEPTH = 0, PVALUE = 0))
  
  keep.mks <- which(mydata$.__enclos_env__$private$pos %in% tot_mks$pos)
  
  pos <- mydata$.__enclos_env__$private$pos[keep.mks]
  depth_Ref_m <- mydata$.__enclos_env__$private$ref[[1]][,keep.mks]
  depth_Alt_m <- mydata$.__enclos_env__$private$alt[[1]][,keep.mks]
  
  depth_Ref <- list(depth_Ref_m)
  depth_Alt <- list(depth_Alt_m)
  config <- mydata$.__enclos_env__$private$config[[1]][keep.mks]
  
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
  real_phase <- real_phases[which(real_phases$pos%in%pos),][,2]
  
  config[which(config==1)] <- "B3.7"
  config[which(config==2 | config==3)] <- "D1.10"
  config[which(config==4 | config==5)] <- "D2.15"
  
  map_info <- data.frame("mk.name"= mydata$.__enclos_env__$private$SNP_Names[keep.mks],
                         "pos" = pos,
                         "rf" = dist.gus,
                         "type"= config,
                         "est.phases"= phases.gus,
                         "real.phases"= real_phase)
  
  return(map_info)
}
# the errors report include markers with distortion and redundants
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
  pos.error <- sapply(strsplit(rownames(meth.error), split = "_"), "[",2)
  ind.error <- paste0(sapply(strsplit(rownames(meth.error), split = "_"), "[", 3), "_", sapply(strsplit(rownames(meth.error), split = "_"), "[", 4))
  meth.error <- as.data.frame(cbind(ind.error, pos.error, meth.error))
  colnames(meth.error) <- c("MK", "POS", "A", "AB", "BA", "B")
  error.info <- merge(gab.geno, meth.geno)
  error.info <- merge(error.info, meth.error)
  return (error.info)
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
