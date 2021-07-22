# Functions
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
  if (is(x$data.name, "outcross")) {
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
    # If there are multiallelic markers, it returns NULL
    multi <- which(sapply(phases, is.null))
    phases[multi] <- NA
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

create_filters_report <- function(onemap_obj, SNPCall, CountsFrom, GenoCall, seed, depth) {
  # onemap_prob <- filter_prob(onemap_obj, threshold = 0.8)
  onemap_mis <- filter_missing(onemap_obj, threshold = 0.25)
  bins <- find_bins(onemap_mis)
  onemap_bins <- create_data_bins(onemap_mis, bins)
  segr <- test_segregation(onemap_bins)
  distorted <- select_segreg(segr, distorted = T)
  no_distorted <- select_segreg(segr, distorted = F, numbers = T)
  twopts <- rf_2pts(onemap_bins) # redundant markers are removed
  seq1 <- make_seq(twopts, no_distorted)
  total_variants <- onemap_obj[[3]]
  filters_tab <- data.frame("higher than 25% missing" = onemap_obj$n.mar - onemap_mis$n.mar,
                            "n_markers"= total_variants,
                            "distorted_markers"=length(distorted), # After red removed
                            "redundant_markers"=total_variants - length(bins[[1]]),
                            "after_filters"= length(seq1$seq.num),
                            "SNPCall" = SNPCall,
                            "GenoCall" = GenoCall,
                            "CountsFrom" = CountsFrom, 
                            seed, 
                            depth)
  
  write_report(filters_tab, "filters_report")
  return(seq1)
}

create_maps_report <- function(input.seq, 
                               tot_mks,
                               gab, 
                               SNPCall, 
                               GenoCall, 
                               fake, 
                               CountsFrom,
                               real_phases, seed, depth, max_cores) {
  
  if(!fake){
    true_mks <- input.seq$seq.num[which(input.seq$data.name$POS[input.seq$seq.num] %in% tot_mks[,2])]
    seq_true <- make_seq(input.seq$twopt, true_mks) # only true markers are mapped
  } else {
    real.mks <- input.seq$data.name$POS[input.seq$seq.num] %in% tot_mks[,2]
    real.mks[which(real.mks == T)] <- "true marker"
    real.mks[which(real.mks == F)] <- "false positive"
    seq_true <- input.seq
  }
  
  if(length(seq_true$seq.num) > 100){
    # batch size close to 60 and the overlap is 3/5 of the size (according with paper)
    #div <- round((length(input.seq$seq.num)/60),0)
    #size = round(length(input.seq$seq.num)/div,0)
    #overlap = round(size*(3/5),0)
    #around = 10
    
    batch_size <- pick_batch_sizes(seq_true,
                                   size = 50,
                                   overlap = 30,
                                   around = 10)
    
    map_df <- map_avoid_unlinked(input.seq = seq_true, size = batch_size, 
                                 phase_cores = max_cores, overlap = 30)
    
  } else {
    map_df <- map_avoid_unlinked(seq_true)
  }
  
  phases <- phaseToOPGP_OM(x = map_df)
  types <- input.seq$data.name$segr.type[map_df[[1]]]
  pos <- input.seq$data.name$POS[map_df[[1]]]
  
  if(!fake){
    real_type <- rep(NA, length(types))
    temp_type <- gab$segr.type[which(as.character(gab$POS) %in% input.seq$data.name$POS[map_df[[1]]])]
    real_type[which(input.seq$data.name$POS[map_df[[1]]] %in% as.character(gab$POS))] <- temp_type
    real_type[which(is.na(real_type))] <- "non-informative"
    real_phase <- real_phases[which(real_phases[,1] %in% input.seq$data.name$POS[map_df[[1]]]),2]
    poscM <- tot_mks$pos.map[which(as.numeric(as.character(tot_mks$pos)) %in% as.numeric(as.character(pos)))]
    poscM.norm <- c(0,cumsum(diff(poscM)))
    diff= sqrt((poscM.norm - c(0,cumsum(haldane(map_df$seq.rf))))^2)
    
    map_info <- data.frame(seed,
                           depth,
                           "mk.name"= colnames(input.seq$data.name$geno)[map_df[[1]]],
                           "pos" = input.seq$data.name$POS[map_df[[1]]],
                           "rf" = c(0,cumsum(haldane(map_df[[3]]))),
                           "type"= types,
                           "real.type" = real_type,
                           "est.phases"= unlist(phases),
                           "real.phases"= real_phase,
                           "real.mks" = "true marker",
                           "SNPCall" = SNPCall,
                           "GenoCall" = GenoCall,
                           "CountsFrom" = CountsFrom,
                           "fake" = "without-false",
                           "poscM" = poscM,
                           "poscM.norm" = poscM.norm,
                           "diff" = diff)
  } else { # Including fake markers is not possible to do the comparisions
    # The fake markers can also be multiallelic markers
    real.mks <- real.mks[input.seq$seq.num %in% map_df$seq.num]
    map_info <- data.frame(seed,
                           depth,
                           "mk.name"= colnames(input.seq$data.name$geno)[map_df[[1]]],
                           "pos" = input.seq$data.name$POS[map_df[[1]]],
                           "rf" = c(0,cumsum(haldane(map_df[[3]]))),
                           "type"= types,
                           "real.type" = NA,
                           "est.phases"= unlist(phases),
                           "real.phases"= NA,
                           "real.mks" = real.mks,
                           "SNPCall" = SNPCall,
                           "GenoCall" = GenoCall,
                           "CountsFrom" = CountsFrom,
                           "fake" = "with-false",
                           "poscM" = NA,
                           "poscM.norm" = NA, 
                           "diff" = NA)
  }
  
  return(list(map_df, map_info))
}


create_gusmap_report <- function(vcf_file, gab, SNPCall, GenoCall, fake, CountsFrom, tot_mks, real_phases, seed, depth){
  ## Maps with gusmap
  RAfile <- VCFtoRA(vcf_file, makePed = T)
  filelist = list.files(pattern = ".*_ped.csv")
  
  ped.file <- read.csv(filelist)
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
  filelist = list.files(pattern = ".*.ra.tab")
  RAdata <- readRA(filelist, pedfile = "ped.file.csv", 
                   filter = list(MAF=0.05, MISS=0.25, BIN=0, DEPTH=0, PVALUE=0.05), sampthres = 0)
  
  mydata <- makeFS(RAobj = RAdata, pedfile = "ped.file.csv", 
                   filter = list(MAF = 0.05, MISS = 0.25,
                                 BIN = 1, DEPTH = 0, PVALUE = 0.05, MAXDEPTH=1000))
  
  # Suggested in vignette
  #mydata$rf_2pt(nClust=1)
  #mydata$createLG()
  #mydata$computeMap(mapped = FALSE, inferOPGP = F) Error
  
  # Alternative way
  if(!fake){
    keep.mks <- which(mydata$.__enclos_env__$private$pos %in% tot_mks$pos)
    pos <- mydata$.__enclos_env__$private$pos[keep.mks]
    depth_Ref_m <- mydata$.__enclos_env__$private$ref[[1]][,keep.mks]
    depth_Alt_m <- mydata$.__enclos_env__$private$alt[[1]][,keep.mks]
    config <- mydata$.__enclos_env__$private$config[[1]]
    # If there are inferred config
    infer <- which(is.na(config))
    if(length(infer) > 0)
      config[infer] <- mydata$.__enclos_env__$private$config_infer[[1]][infer]
    config <- config[keep.mks]
  } else {
    pos <- mydata$.__enclos_env__$private$pos
    depth_Ref_m <- mydata$.__enclos_env__$private$ref[[1]]
    depth_Alt_m <- mydata$.__enclos_env__$private$alt[[1]]
    config <- mydata$.__enclos_env__$private$config[[1]]
    # If there are inferred config
    infer <- which(is.na(config))
    if(length(infer) > 0)
      config[infer] <- mydata$.__enclos_env__$private$config_infer[[1]][infer]
    real.mks <- mydata$.__enclos_env__$private$pos %in% tot_mks$pos
    real.mks[which(real.mks == T)] <- "true marker"
    real.mks[which(real.mks == F)] <- "false positive"
  }
  
  # Automatically remove markers with pos = inf
  inf.mk <- 0
  rast.pos <- pos
  while(length(inf.mk) > 0){
    depth_Ref <- list(depth_Ref_m)
    depth_Alt <- list(depth_Alt_m)
    
    phases.gus <- GUSMap:::infer_OPGP_FS(depth_Ref_m, depth_Alt_m, 
                                         config, ep=0.001, reltol=1e-3)
    
    rf_est <- GUSMap:::rf_est_FS(init_r = 0.01, ep = 0.001, 
                                 ref = depth_Ref, 
                                 alt = depth_Alt, 
                                 OPGP=list(as.integer(phases.gus)),
                                 nThreads = 1)
    
    # Remove Inf markers
    inf.mk <- which(is.infinite(haldane(rf_est$rf)))
    if(length(inf.mk) > 0){
      depth_Ref_m <- depth_Ref_m[,-inf.mk]
      depth_Alt_m <- depth_Alt_m[,-inf.mk]
      config <- config[-inf.mk]
      rast.pos <- rast.pos[-inf.mk]
    }
  }
  
  dist.gus <- c(0,cumsum(haldane(rf_est$rf))) # haldane mapping function - 
  # I checked with gusmap example 
  # doing as suggested in vignette
  phases.gus[which(phases.gus == 1 | phases.gus == 4)] <- 17
  phases.gus[which(phases.gus == 2 | phases.gus == 3)] <- 18
  phases.gus[which(phases.gus == 5 | phases.gus == 8)] <- 19
  phases.gus[which(phases.gus == 6 | phases.gus == 7)] <- 20
  phases.gus[which(phases.gus == 9 | phases.gus == 12)] <- 21
  phases.gus[which(phases.gus == 10 | phases.gus == 11)] <- 22
  phases.gus[which(phases.gus == 13 | phases.gus == 16)] <- 23
  phases.gus[which(phases.gus == 14 | phases.gus == 15)] <- 24
  real_phase <- real_phases[which(real_phases$pos%in%rast.pos),][,2]
  
  config[which(config==1)] <- "B3.7"
  config[which(config==2 | config==3)] <- "D1.10"
  config[which(config==4 | config==5)] <- "D2.15"
  
  real_type <- rep(NA, length(config))
  temp_type <- gab$segr.type[which(gab$POS %in% rast.pos)]
  real_type[which(rast.pos %in% as.character(gab$POS))] <- temp_type
  real_type[which(is.na(real_type))] <- "non-informative"
  poscM <- tot_mks$pos.map[which(as.numeric(as.character(tot_mks$pos)) %in% as.numeric(as.character(rast.pos)))]
  poscM.norm <- poscM-poscM[1]
  diff= sqrt((poscM.norm - dist.gus)^2)
  
  if(!fake){
    map_info <- data.frame(seed,
                           depth,
                           "mk.name"= mydata$.__enclos_env__$private$SNP_Names[keep.mks][which(pos%in%rast.pos)],
                           "pos" = rast.pos,
                           "rf" = dist.gus,
                           "type"= config,
                           "real.type" = real_type,
                           "est.phases"= phases.gus,
                           "real.phases"= real_phase,
                           "real.mks" = "true marker",
                           "SNPCall" = SNPCall,
                           "GenoCall" = GenoCall,
                           "CountsFrom" = CountsFrom,
                           "fake" = fake,
                           "poscM" = poscM,
                           "poscM.norm" = poscM.norm,
                           "diff" = diff)
  } else {
    map_info <- data.frame(seed,
                           depth,
                           "mk.name"= mydata$.__enclos_env__$private$SNP_Names[which(pos%in%rast.pos)],
                           "pos" = rast.pos,
                           "rf" = dist.gus,
                           "type"= config,
                           "real.type" = NA,
                           "est.phases"= phases.gus,
                           "real.phases"= NA,
                           "real.mks" = real.mks[which(pos%in%rast.pos)],
                           "SNPCall" = SNPCall,
                           "GenoCall" = GenoCall,
                           "CountsFrom" = CountsFrom,
                           "fake" = fake,
                           "poscM" = NA,
                           "poscM.norm" = NA,
                           "diff" = NA)
  }
  
  map_df <- mydata
  return(list(map_df, map_info))
}


# deprecated
# the errors report include markers with distortion and redundants
create_errors_report <- function(onemap_obj, gab, SNPCall, GenoCall, CountsFrom, seed, depth) {
  pos <- which(gab[[9]] %in% onemap_obj[[9]])
  if(length(pos) < 1) {
    out_data <- as.data.frame(matrix(NA, nrow=2, ncol=11))
  } else {
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
    out_data <- data.frame(seed, depth, SNPCall, GenoCall, CountsFrom, error.info)
  }
  outname <- paste0("errors_report")
  write_report(out_data, outname)
}

write_report <- function(filters_tab, out_name, max_cores=1) {
  vroom::vroom_write(filters_tab, paste0(out_name, ".tsv.gz"), num_threads= max_cores)
}

# Deprecated
make_vcf <- function(vcf.old, depths, allele_file, out_vcf, cores=3){
  # The input od polyRAD need to be a VCF, then this part takes the allele depth from "depths" and put at AD field of input vcf
  idx <- system(paste0("grep -in 'CHROM' ", vcf.old), intern = T) # This part only works in linux OS
  idx.i <- strsplit(idx, split = ":")[[1]][1]
  seed <- sample(1:10000, 1)
  system(paste0("head -n ", idx.i," ", vcf.old, " > head.",seed))
  
  vcf.tab <- read.table(vcf.old, stringsAsFactors = F)
  
  if(all(rownames(depths[[1]]) %in% paste0(vcf.tab[,1], "_", vcf.tab[,2]))){
    
    vcf.init <- vcf.tab[,1:8]
    AD.colum <- rep("GT:AD", dim(vcf.init)[1])
    vcf.init <- cbind(vcf.init, AD.colum)
    
    rs <- rownames(depths[[1]])
    vcf.init[,3] <- rs
  } else {
    temp.tab <- read.table(allele_file)
    vcf.init <- cbind(temp.tab[,1:2],
                      paste0(temp.tab[,1], "_", temp.tab[,2]), 
                      temp.tab[,3:4],
                      rep(".", dim(temp.tab)[1]),
                      rep(".", dim(temp.tab)[1]), 
                      rep(".", dim(temp.tab)[1]), 
                      rep("GT:AD",dim(temp.tab)[1]))
  }
  
  ind.n <- colnames(depths[[1]]) # The names came in different order
  
  header <- strsplit(idx, split = "\t")[[1]]
  ind.vcf <- header[10:length(header)]
  ind.n <- factor(ind.n, levels = ind.vcf)
  
  depths[[1]] <- depths[[1]][,order(ind.n)]
  depths[[2]] <- depths[[2]][,order(ind.n)]
  
  osize <- as.matrix(depths[[1]] + depths[[2]])
  # Remove missing data
  rm.mk <- which(apply(osize, 1, function(x) all(x==0)))
  if(length(rm.mk) > 0){
    depths[[1]] <- depths[[1]][-rm.mk,]
    depths[[2]] <- depths[[2]][-rm.mk,]
    vcf.init <- vcf.init[-rm.mk,]
  }
  osize <- as.matrix(depths[[1]] + depths[[2]])
  oref <- as.matrix(depths[[1]])
  
  # Getting genotypes using updog normal model
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl = cl)
  gene_est <- foreach(i = 1:nrow(osize)) %dopar% {
    ## fit flexdog
    fout <- updog::flexdog(refvec  = oref[i,],
                           sizevec = osize[i,],
                           ploidy  = 2,
                           model = "norm")
    fout
  }
  
  parallel::stopCluster(cl)
  geno <- t(sapply(gene_est, "[[", 10))
  geno[which(geno == 1)] <- "0/1"
  geno[which(geno == 0)] <- "1/1"
  geno[which(geno == 2)] <- "0/0" # Updog counts by the reference allele
  
  comb.depth <- matrix(paste0(geno, ":", as.matrix(depths[[1]]), ",", as.matrix(depths[[2]])), ncol = ncol(depths[[2]]))
  colnames(comb.depth) <- ind.vcf
  #hmc.file <- cbind(rs, comb.depth)
  
  vcf.body <- cbind(vcf.init, comb.depth)
  
  write.table(vcf.body, file = paste0("temp.body.", seed), quote = FALSE, sep = "\t", row.names = FALSE, col.names = F)
  
  system(paste0("cat head.",seed," temp.body.",seed," > ", out_vcf))
  return(out_vcf)
}

# deprecated
adapt2app <- function(datas_up_inp){
  
  # For errors
  datas_up_inp[[1]]$ref[which(datas_up_inp[[1]]$ref == ".")] <- NA
  datas_up_inp[[1]]$alt[which(datas_up_inp[[1]]$alt == ".")] <- NA
  datas_up_inp[[1]]$ref <- as.numeric(datas_up_inp[[1]]$ref)
  datas_up_inp[[1]]$alt <- as.numeric(datas_up_inp[[1]]$alt)
  
  datas_up_inp[[1]] <- fix_genocall_names(datas_up_inp[[1]])
  
  ####
  datas_up_inp[[2]] <- datas_up_inp[[2]][,-3]
  colnames(datas_up_inp[[2]]) <- c("seed", "depth", "mk.name", "pos" , "rf" , "type", "real.type", 
                                   "est.phases", "real.phases", "real.mks", "SNPCall", "GenoCall", 
                                   "CountsFrom", "fake", "poscM", "poscM.norm", "diff")
  
  datas_up_inp[[2]]$fake[which(datas_up_inp[[2]]$real.mks ==99)] <- "without-false"
  datas_up_inp[[2]]$fake[which(datas_up_inp[[2]]$real.mks < 2)] <- "with-false"
  
  datas_up_inp[[2]]$real.mks[which(datas_up_inp[[2]]$real.mks == 99 | datas_up_inp[[2]]$real.mks == 1)] <- "true markers"
  datas_up_inp[[2]]$real.mks[which(datas_up_inp[[2]]$real.mks == 0)] <- "false positives"
  
  datas_up_inp[[2]] <- fix_genocall_names(datas_up_inp[[2]])
  
  ###
  colnames(datas_up_inp[[3]]) <- c("seed", "depth", "mis_markers", "n_markers", "distorted_markers", "redundant_markers", "SNPCall", "GenoCall", "CountsFrom")
  
  datas_up_inp[[3]] <- fix_genocall_names(datas_up_inp[[3]])
  
  ###
  datas_up_inp[[4]] <- as.data.frame(datas_up_inp[[4]])
  datas_up_inp[[4]]$fake <- as.character(datas_up_inp[[4]]$fake)
  datas_up_inp[[4]]$fake[which(datas_up_inp[[4]]$fake =="FALSE")] <- "without-false"
  datas_up_inp[[4]]$fake[which(datas_up_inp[[4]]$fake =="TRUE")] <- "with-false"
  datas_up_inp[[4]]$fake <- as.factor(datas_up_inp[[4]]$fake)
  
  datas_up_inp[[4]] <- fix_genocall_names(datas_up_inp[[4]])
  
  #####
  datas_up_inp[[5]] <- datas_up_inp[[5]][,-7] ## Ajustar
  names(datas_up_inp[[5]]) <- c("depth", "seed", "SNPCall", "(1)", "(2)", "(3)", "(4)", "(5)") ## Ajustar
  datas_up_inp[[5]] <- gather(datas_up_inp[[5]], key,value,-SNPCall, -depth, -seed)
  
  datas_up_inp[[3]]$GenoCall <- factor(datas_up_inp[[3]]$GenoCall)
  datas_up_inp[[3]] <- gather(datas_up_inp[[3]],key,value, -CountsFrom, -seed, -depth, -SNPCall, -GenoCall)
  datas_up_inp[[3]]$key <- factor(datas_up_inp[[3]]$key)
  
  # Defining the options
  cout <- table(datas_up_inp[[1]]$seed, datas_up_inp[[1]]$depth)
  depthNames <- colnames(cout)
  depths <- seeds <- seedsNames <- vector()
  for(i in 1:length(depthNames)){
    for(j in 1:nrow(cout)){
      if(cout[j,i] > 0){
        temp <- rownames(cout)[j]
        seedsNames <- c(seedsNames, paste0("Depth ", depthNames[i], " seed ", temp))
        depths <- c(depths, depthNames[i])
        seeds <- c(seeds, temp)
      }
    }
  }
  
  depth_choice <- as.list(unique(datas_up_inp[[1]]$depth))
  names(depth_choice) <- as.character(unique(datas_up_inp[[1]]$depth))
  seeds_choice <- as.list(1:length(seedsNames))
  names(seeds_choice) <- as.character(seedsNames)
  
  result_list <- list("data1" = datas_up_inp[[1]], "data2"= datas_up_inp[[2]], 
                      "data3"=datas_up_inp[[3]], "data4"=datas_up_inp[[4]], 
                      "data5"=datas_up_inp[[5]], "choices" = list(depths, seeds, seeds_choice, depth_choice))
  
  return(result_list)
}

#deprecated
fix_genocall_names <- function(data_broken){
  data_broken$GenoCall <- as.factor(data_broken$GenoCall)
  temp <- levels(data_broken$GenoCall)
  idx <- match(c("default", "default0.05"), temp)
  temp[idx] <- c("OneMap_version2", "SNPCaller0.05")
  data_broken$GenoCall <- as.character(data_broken$GenoCall)
  data_broken$GenoCall[data_broken$GenoCall == "default"] <- "OneMap_version2"
  data_broken$GenoCall[data_broken$GenoCall == "default0.05"] <- "SNPCaller0.05"
  data_broken$GenoCall <- factor(data_broken$GenoCall, labels = temp, levels = temp)
  return(data_broken)
}

update_fake_info <- function(info_fake, simu_onemap_obj, ref_alt_alleles, simulated_phases){
  info_correct <- info_fake
  est.pos <- info_fake[[2]]$pos
  real.type <- rep(NA, nrow(info_correct[[2]]))
  temp.type <- simu_onemap_obj$segr.type[which(simu_onemap_obj$POS %in% est.pos)]
  real.type[which(est.pos %in% as.character(simu_onemap_obj$POS))] <- temp.type
  real.type[which(is.na(real.type))] <- "non-informative"
  poscM <- ref_alt_alleles$pos.map[which(as.numeric(as.character(ref_alt_alleles$pos)) %in% as.numeric(as.character(est.pos)))]
  poscM.norm <- poscM-poscM[1]
  diff <- sqrt((poscM.norm - info_fake[[2]]$rf)^2)
  real.phase <- simulated_phases[which(simulated_phases$pos%in%est.pos),][,2]
  
  info_correct[[2]]$real.type <- real.type
  info_correct[[2]]$real.phases <- real.phase
  info_correct[[2]]$fake <- FALSE
  info_correct[[2]]$poscM <- poscM
  info_correct[[2]]$poscM.norm <- poscM.norm
  info_correct[[2]]$diff <- diff
  
  return(info_correct)
}
