# Functions
create_map_report <- function(input.seq, CountsFrom, SNPCall, GenoCall, max_cores){
  # Check genome position
  pos <- as.numeric(input.seq$data.name$POS[input.seq$seq.num])
  sort.pos <- sort(pos)
  if(all(pos != sort.pos)){
    cat("The markers are not ordered by genome position")
    input.seq <- make_seq(input.seq$twopt, input.seq$seq.num[order(as.numeric(input.seq$data.name$POS[input.seq$seq.num]))])
  } 
  
  if(length(input.seq$seq.num) > 100){
    # batch size close to 60 and the overlap is 3/5 of the size (according with paper)
    #div <- round((length(input.seq$seq.num)/60),0)
    #size = round(length(input.seq$seq.num)/div,0)
    #overlap = round(size*(3/5),0)
    #around = 10
    
    batch_size <- pick_batch_sizes(input.seq,
                                   size = 50,
                                   overlap = 30,
                                   around = 10)
    
    map_out <- map_avoid_unlinked(input.seq, 
                                  size = batch_size, 
                                  phase_cores = max_cores, 
                                  overlap = 30)
  } else {
    map_out <- map_avoid_unlinked(input.seq)
  }
  
  sizes_df <- data.frame(CountsFrom, SNPCall, GenoCall, "mks" = colnames(map_out$data.name$geno)[map_out$seq.num],
                         "pos" = map_out$data.name$POS[map_out$seq.num], rf = cumsum(c(0,kosambi(map_out$seq.rf))),
                         type = map_out$data.name$segr.type[map_out$seq.num], phases = phaseToOPGP_OM(map_out))
  return(list(map_out, sizes_df))
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
  if (is(x$data.name ,"outcross")) {
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

create_filters_report <- function(onemap_obj, SNPCall,CountsFrom, GenoCall, chromosome) {
  # onemap_prob <- filter_prob(onemap_obj, threshold = 0.8)
  onemap_mis <- filter_missing(onemap_obj, threshold = 0.25)
  onemap_mis <- onemap::filter_missing(onemap_obj, threshold = 0.25)
  bins <- onemap::find_bins(onemap_mis)
  onemap_bins <- create_data_bins(onemap_mis, bins)
  twopts <- rf_2pts(input.obj = onemap_bins, rm_mks = T) # Do not keep redundant markers
  new_obj <- twopts$data.name
  segr <- onemap::test_segregation(new_obj)
  distorted <- onemap::select_segreg(segr, distorted = T)
  no_distorted <- onemap::select_segreg(segr, distorted = F, numbers = T)
  chr <- which(new_obj$CHROM %in% chromosome)
  chr_no_dist <- chr[which(chr%in%no_distorted)]
  seq1 <- make_seq(twopts, chr_no_dist)
  total_variants <- onemap_obj[[3]]
  lgs <- group(seq1)
  if(all(lgs$groups == 0)) {
    nongroup <- length(seq1$seq.num)
    lg1 <- seq1
  } else {
    lg1 <- make_seq(lgs, as.numeric(names(which.max(table(lgs$groups[-which(lgs$groups== 0)])))))
    nongroup <- length(seq1$seq.num) - length(lg1$seq.num)
  }
  filters_tab <- data.frame(CountsFrom,
                            SNPCall,
                            GenoCall,
                            "higher than 25% missing" = onemap_obj$n.mar - onemap_mis$n.mar,
                            "n_markers"= total_variants,
                            "n_markers_selected_chr" = length(chr),
                            "selected_chr_no_dist" = length(chr_no_dist),
                            "distorted_markers"= length(distorted),
                            "redundant_markers"= total_variants - length(bins[[1]]),
                            "non-grouped_markers" = nongroup)
  write_report(filters_tab, paste0(SNPCall, "_",GenoCall, "_",CountsFrom,"_filters_report"))
  return(lg1)
}

write_report <- function(tab, out_name, max_cores=1) {
  vroom::vroom_write(tab, paste0(out_name, ".tsv.gz"), num_threads= max_cores)
}

make_vcf <- function(vcf.old, depths, method, allele_file, out_vcf){
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
    temp.tab <- read.table(allele_file)
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
  
  system(paste0("cat head.",seed," temp.body.",seed," > ", out_vcf))
  return(out_vcf)
}

create_gusmap_report <- function(vcf_file,SNPCall, CountsFrom, GenoCall, parent1, parent2){
  ## Maps with gusmap
  RAfile <- VCFtoRA(vcf_file, makePed = T)
  filelist = list.files(pattern = ".*_ped.csv")
  
  ped.file <- read.csv(filelist)
  ID <- c(1:(dim(ped.file)[1]))
  idx.P1 <- which(as.character(ped.file$SampleID) == parent1)
  idx.P2 <- which(as.character(ped.file$SampleID) == parent2)
  
  mother <- c(rep(idx.P1,dim(ped.file)[1]))
  father <- c(rep(idx.P2,dim(ped.file)[1]))
  fam <- c(rep("F1", dim(ped.file)[1]))
  idx.P1 <- which(as.character(ped.file$SampleID) == parent1)
  idx.P2 <- which(as.character(ped.file$SampleID) == parent2)
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
  
  pos <- mydata$.__enclos_env__$private$pos
  depth_Ref_m <- mydata$.__enclos_env__$private$ref[[1]]
  depth_Alt_m <- mydata$.__enclos_env__$private$alt[[1]]
  config <- mydata$.__enclos_env__$private$config[[1]]
  # If there are inferred config
  infer <- which(is.na(config))
  if(length(infer) > 0)
    config[infer] <- mydata$.__enclos_env__$private$config_infer[[1]][infer]
  
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
  
  #rf_est$rf[which(rf_est$rf > 0.5)] <- 0.4999999
  dist.gus <- c(0,cumsum(haldane(rf_est$rf)))
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
  
  file.name <- paste0(SNPCall, "_", CountsFrom, "_", GenoCall)
  map_out <- mydata
  
  map_info <- data.frame("CountsFrom" = CountsFrom,
                         "SNPCall"= SNPCall,
                         "GenoCall" = GenoCall,
                         "mks"= mydata$.__enclos_env__$private$SNP_Names[which(pos%in%rast.pos)],
                         "pos" = rast.pos,
                         "rf" = dist.gus,
                         "type"= config,
                         "phases"= phases.gus)
  
  return(list(map_out, map_info))
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
