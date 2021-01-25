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

create_filters_report <- function(onemap_obj, SNPcall, CountsFrom, Genocall) {
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
                            "SNPcall" = SNPcall,
                            "Genocall" = Genocall,
                            "CountsFrom" = CountsFrom)
  
  write_report(filters_tab, paste0("filters_", SNPcall, "_", CountsFrom, "_",Genocall, ".txt"))
  return(seq1)
}

create_maps_report <- function(input.seq, 
                               tot_mks,
                               gab, 
                               SNPcall, 
                               Genocall, 
                               fake, 
                               CountsFrom,
                               real_phases) {
  
  if(!fake){
    true_mks <- input.seq$seq.num[which(input.seq$data.name$POS[input.seq$seq.num] %in% tot_mks[,2])]
    seq_true <- make_seq(input.seq$twopt, true_mks) # only true markers are mapped
  } else {
    real.mks <- input.seq$data.name$POS[input.seq$seq.num] %in% tot_mks[,2]
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
                                 phase_cores = 4, overlap = overlap)
    
  } else {
    map_df <- map_avoid_unlinked(seq_true)
  }
  
  phases <- phaseToOPGP_OM(x = map_df)
  types <- input.seq$data.name$segr.type[map_df[[1]]]
  real_type <- rep(NA, length(types))
  temp_type <- gab$segr.type[which(as.character(gab$POS) %in% input.seq$data.name$POS[map_df[[1]]])]
  real_type[which(input.seq$data.name$POS[map_df[[1]]] %in% as.character(gab$POS))] <- temp_type
  real_type[which(is.na(real_type))] <- "non-informative"
  real_phase <- real_phases[which(real_phases[,1] %in% input.seq$data.name$POS[map_df[[1]]]),2]
  pos <- input.seq$data.name$POS[map_df[[1]]]
  poscM <- tot_mks$pos.map[which(as.numeric(as.character(tot_mks$pos)) %in% as.numeric(as.character(pos)))]
  poscM.norm <- poscM-poscM[1]
  diff= sqrt((poscM.norm - c(0,cumsum(haldane(map_df[[3]]))))^2)
  
  if(!fake){
    map_info <- data.frame("mk.name"= colnames(input.seq$data.name$geno)[map_df[[1]]],
                           "pos" = input.seq$data.name$POS[map_df[[1]]],
                           "rf" = c(0,cumsum(haldane(map_df[[3]]))),
                           "type"= types,
                           "real.type" = real_type,
                           "est.phases"= unlist(phases),
                           "real.phases"= real_phase,
                           "real.mks" = 99,
                           "SNPcall" = SNPcall,
                           "Genocall" = Genocall,
                           "CountsFrom" = CountsFrom,
                           "fake" = fake,
                           "poscM" = poscM,
                           "poscM.norm" = poscM.norm,
                           "diff" = diff)
  } else { # Including fake markers is not possible to do the comparisions
    # The fake markers can also be multiallelic markers
    real.mks <- real.mks[input.seq$seq.num %in% map_df$seq.num]
    map_info <- data.frame("mk.name"= colnames(input.seq$data.name$geno)[map_df[[1]]],
                           "pos" = input.seq$data.name$POS[map_df[[1]]],
                           "rf" = c(0,cumsum(haldane(map_df[[3]]))),
                           "type"= types,
                           "real.type" = NA,
                           "est.phases"= unlist(phases),
                           "real.phases"= NA,
                           "real.mks" = real.mks,
                           "SNPcall" = SNPcall,
                           "Genocall" = Genocall,
                           "CountsFrom" = CountsFrom,
                           "fake" = fake,
                           "poscM" = NA,
                           "poscM.norm" = NA, 
                           "diff" = NA)
  }
  
  save(map_df, file= paste0("map_", SNPcall, "_", CountsFrom, "_",Genocall, "_", fake, ".RData"))
  write_report(map_info, paste0("map_", SNPcall, "_", CountsFrom, "_",Genocall, "_",fake, ".txt"))
  return(map_info)
}


create_gusmap_report <- function(vcf_file, gab, SNPcall, Genocall, fake, CountsFrom, tot_mks, real_phases){
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
                   filter = list(MAF=0.05, MISS=0.5, BIN=0, DEPTH=0, PVALUE=0.05), sampthres = 0)
  
  mydata <- makeFS(RAobj = RAdata, pedfile = "ped.file.csv", 
                   filter = list(MAF = 0.05, MISS = 0.5,
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
    map_info <- data.frame("mk.name"= mydata$.__enclos_env__$private$SNP_Names[keep.mks][which(pos%in%rast.pos)],
                           "pos" = rast.pos,
                           "rf" = dist.gus,
                           "type"= config,
                           "real.type" = real_type,
                           "est.phases"= phases.gus,
                           "real.phases"= real_phase,
                           "real.mks" = 99,
                           "SNPcall" = SNPcall,
                           "Genocall" = Genocall,
                           "CountsFrom" = CountsFrom,
                           "fake" = fake,
                           "poscM" = poscM,
                           "poscM.norm" = poscM.norm,
                           "diff" = diff)
  } else {
    map_info <- data.frame("mk.name"= mydata$.__enclos_env__$private$SNP_Names[which(pos%in%rast.pos)],
                           "pos" = rast.pos,
                           "rf" = dist.gus,
                           "type"= config,
                           "real.type" = NA,
                           "est.phases"= phases.gus,
                           "real.phases"= NA,
                           "real.mks" = real.mks[which(pos%in%rast.pos)],
                           "SNPcall" = SNPcall,
                           "Genocall" = Genocall,
                           "CountsFrom" = CountsFrom,
                           "fake" = fake,
                           "poscM" = NA,
                           "poscM.norm" = NA,
                           "diff" = NA)
  }
  
  outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
  map_df <- mydata
  save(map_df, file = paste0(outname,".RData"))
  write_report(map_info, paste0(outname, ".txt"))
}
# the errors report include markers with distortion and redundants
create_errors_report <- function(onemap_obj, gab, SNPcall, Genocall, CountsFrom) {
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
    out_data <- data.frame(SNPcall, Genocall, CountsFrom, error.info)
  }
  outname <- paste0("errors_", SNPcall, "_", CountsFrom, "_", Genocall)
  write_report(out_data, paste0(outname,".txt"))
}

write_report <- function(filters_tab, out_name) {
  write.table(filters_tab, file=out_name, row.names=F, quote=F, col.names = F)
}


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

adapt2app <- function(data){
  
  # For errors
  data[[1]]$errors <- apply(data[[1]][,10:13], 1, function(x) if(all(x==1)) NA else 1 - x[which.max(x)])
  data[[1]]$ref[which(data[[1]]$ref == ".")] <- NA
  data[[1]]$ref <- as.numeric(data[[1]]$ref)
  # The onemap genotype codification do not diferenciate the homozyotes for each parent
  data[[1]]$gabGT[data[[1]]$gabGT == 3 | data[[1]]$gabGT == 1] <- "homozygous"
  data[[1]]$gabGT[data[[1]]$gabGT == 2] <- "heterozygote"
  data[[1]]$methGT[data[[1]]$methGT == "3" | data[[1]]$methGT == "1"] <- "homozygous"
  data[[1]]$methGT[data[[1]]$methGT == "2"] <- "heterozygote"
  data[[1]]$methGT[data[[1]]$methGT == "0"] <- "missing"
  
  data[[1]]$methGT <- factor(data[[1]]$methGT, labels = c("missing", "homozygous", "heterozygote"), levels = c("missing", "homozygous", "heterozygote"))
  data[[1]]$gabGT <- factor(data[[1]]$gabGT, labels = c("missing", "homozygous", "heterozygote"), levels = c("missing", "homozygous", "heterozygote"))
  
  data[[1]] <- fix_genocall_names(data[[1]])
  
  ####
  data[[2]] <- data[[2]][,-3]
  colnames(data[[2]]) <- c("seed", "depth", "mk.name", "pos" , "rf" , "type", "real.type", 
                           "est.phases", "real.phases", "real.mks", "SNPCall", "GenoCall", 
                           "CountsFrom", "fake", "poscM", "poscM.norm", "diff")
  
  data[[2]]$fake[which(data[[2]]$real.mks ==99)] <- "without-false"
  data[[2]]$fake[which(data[[2]]$real.mks < 2)] <- "with-false"
  
  data[[2]]$real.mks[which(data[[2]]$real.mks == 99 | data[[2]]$real.mks == 1)] <- "true markers"
  data[[2]]$real.mks[which(data[[2]]$real.mks == 0)] <- "false positives"
  
  data[[2]] <- fix_genocall_names(data[[2]])
  
  ###
  colnames(data[[3]]) <- c("seed", "depth", "mis_markers", "n_markers", "distorted_markers", "redundant_markers", "SNPCall", "GenoCall", "CountsFrom")
  
  data[[3]] <- fix_genocall_names(data[[3]])
  
  ###
  data[[4]] <- as.data.frame(data[[4]])
  data[[4]]$fake <- as.character(data[[4]]$fake)
  data[[4]]$fake[which(data[[4]]$fake =="FALSE")] <- "without-false"
  data[[4]]$fake[which(data[[4]]$fake =="TRUE")] <- "with-false"
  data[[4]]$fake <- as.factor(data[[4]]$fake)
  
  data[[4]] <- fix_genocall_names(data[[4]])
  
  #####
  data[[5]] <- data[[5]][,-7] ## Ajustar
  names(data[[5]]) <- c("depth", "seed", "SNPCall", "(1)", "(2)", "(3)", "(4)", "(5)") ## Ajustars
  data[[5]] <- gather(data[[5]], key,value,-SNPCall, -depth, -seed)
  
  data[[3]]$GenoCall <- factor(data[[3]]$GenoCall)
  data[[3]] <- gather(data[[3]],key,value, -CountsFrom, -seed, -depth, -SNPCall, -GenoCall)
  data[[3]]$key <- factor(data[[3]]$key)
  
  # Defining the options
  cout <- table(data[[1]]$seed, data[[1]]$depth)
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
  
  depth_choice <- as.list(unique(data[[1]]$depth))
  names(depth_choice) <- as.character(unique(data[[1]]$depth))
  seeds_choice <- as.list(1:length(seedsNames))
  names(seeds_choice) <- as.character(seedsNames)
  
  result_list <- list("data1" = data[[1]], "data2"= data[[2]], 
                      "data3"=data[[3]], "data4"=data[[4]], 
                      "data5"=data[[5]], "choices" = list(depths, seeds, seeds_choice, depth_choice))
  
  return(result_list)
}

fix_genocall_names <- function(data){
  data$GenoCall <- as.factor(data$GenoCall)
  temp <- levels(data$GenoCall)
  idx <- grep(c("default", "default0.05"), temp)
  temp[idx] <- c("OneMap_version2", "SNPCaller0.05")
  data$GenoCall <- as.character(data$GenoCall)
  data$GenoCall[data$GenoCall == "default"] <- "OneMap_version2"
  data$GenoCall[data$GenoCall == "default0.05"] <- "SNPCaller0.05"
  data$GenoCall <- factor(data$GenoCall, labels = temp, levels = temp)
  return(data)
}

