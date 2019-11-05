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
  segr <- onemap::test_segregation(onemap_obj)
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
