##################################################################################
# Functions to create PedigreeSim inputs from 
# * A unphased VCF file (haplotypes are simulated)
# * Pirs diploid simulated SNPs and indels - needs optimization in WDL
#
# Genetic distances simulated with:
# * With a pre-defined cM by Mb rate
# * with a reference linkage map
# 
# Written by Cristiane Taniguti
##################################################################################

#' Creates PedigreeSim input files from empirical non-phased VCF file
#' The same alleles and proportion of marker types are kept from the VCF and the haplotypes are simulated
#' This function needs adaptation to deal with other populations and ploidy
#' It is important that MNP markers are in splitted format
create_haplo <- function(vcfR.obj, ref.map, seed, P1, P2, filename = "founders.txt"){
  set.seed(seed)
  GT <- vcfR.obj@gt
  inds <- colnames(GT)
  parents <- data.frame("P1" = GT[,which(inds==P1)], "P2" = GT[,which(inds==P2)])
  GT_matrix <- matrix(sapply(strsplit(as.matrix(parents), ":"), "[[",1), ncol = 2)
  
  # Only calculate the percentage of non-informative, D1.10, D2.15 and B3.7 markers and 
  # after the haplotypes are randomly choosed
  # This function do not consider phased genotypes
  if(any(grepl("[|]", GT_matrix))){
    GT_matrix <- gsub("[|]", "/", as.matrix(GT_matrix))
    GT_matrix[which(GT_matrix == "1/0")] <- "0/1"
    GT_matrix[which(GT_matrix == "2/0")] <- "0/2"
    GT_matrix[which(GT_matrix == "3/0")] <- "0/3"
    GT_matrix[which(GT_matrix == "2/1")] <- "1/2"
    GT_matrix[which(GT_matrix == "3/1")] <- "1/3"
    GT_matrix[which(GT_matrix == "3/2")] <- "2/3"
  }
  
  # keep only biallelic
  rm_multi <- which(apply(GT_matrix, 1, function(x) any(grepl("2", x))))
  
  # remove markers out of ref.map interval
  pos.vcf <- as.numeric(as.character(vcfR.obj@fix[,2]))
  ref.alleles <- vcf@fix[,4]
  alt.alleles <- vcf@fix[,5]
  mk.names <- unique(paste0(vcf@fix[,1], "_", vcf@fix[,2]))
  
  rm.mks <- unique(which(pos.vcf < min(ref.map$bp)), which(pos.vcf > max(ref.map$bp)))
  
  rm.mks <- unique(c(rm_multi, rm.mks))
  if(length(rm.mks) > 0){
    GT_matrix <- GT_matrix[-rm.mks,]
    ref.alleles <- ref.alleles[-rm.mks]
    alt.alleles <- alt.alleles[-rm.mks]
    mk.names <- mk.names[-rm.mks]
  }
  n.mk <- nrow(GT_matrix)
  
  mk.type <- mk.type.num <- rep(NA, n.mk)
  
  P1_1 <- sapply(strsplit(GT_matrix[,1], "/"), "[", 1)
  P1_2 <- sapply(strsplit(GT_matrix[,1], "/"), "[", 2)
  P2_1 <- sapply(strsplit(GT_matrix[,2], "/"), "[", 1)
  P2_2 <- sapply(strsplit(GT_matrix[,2], "/"), "[", 2)
  
  # Marker types
  GT_parents <- cbind(P1_1, P1_2,P2_1, P2_2)
  idx <- which(P1_1 == "." | P2_1 == "." |  P1_2 == "." | P2_2 == ".")
  GT_parents[idx,] <- NA
  
  idx <- apply(GT_parents, 1, function(x) (length(x) -2) == length(unique(x)))
  idx.sub <- which(P1_1[idx] != P1_2[idx] & P2_1[idx] != P2_2[idx])
  mk.type[idx][idx.sub] <- "B3.7"
  mk.type.num[idx][idx.sub] <- 4
  idx.sub <- which(P1_1[idx] != P1_2[idx] & P2_1[idx] == P2_2[idx])
  mk.type[idx][idx.sub] <- "D1.10"
  mk.type.num[idx][idx.sub] <- 6
  idx.sub <- which(P1_1[idx] == P1_2[idx] & P2_1[idx] != P2_2[idx])
  mk.type[idx][idx.sub] <- "D2.15"
  mk.type.num[idx][idx.sub] <- 7
  
  # Non-informative markers are also included
  
  idx <- which(GT_parents[,1] == GT_parents[,2] & GT_parents[,3] == GT_parents[,4])
  mk.type[idx] <- "homo-non-info"
  
  # Here we consider that the marker types are equally distributed by the genome
  mk.type.nonna <- mk.type[!is.na(mk.type)]
  
  doses <- round(c(sum(mk.type.nonna == "homo-non-info")/length(mk.type.nonna)*100, 
                   sum(mk.type.nonna == "D1.10")/length(mk.type.nonna)*100, 
                   sum(mk.type.nonna == "B3.7")/length(mk.type.nonna)*100,
                   sum(mk.type.nonna == "D2.15")/length(mk.type.nonna)*100), 0)
  doses <- doses[-1]
  doses <- c(100-sum(doses), doses)
  
  type.quant <- list()
  for(i in 1:(length(doses)-1)){
    type.quant[[i]] <- round((doses[i]/100)*length(mk.names),0)
  }
  type.quant[[i+1]] <- length(mk.names) - sum(unlist(type.quant))
  
  # Defining marker types distribuition
  pos.mks <- 1:length(mk.names)
  pos.non <- sample(pos.mks, type.quant[[1]])
  pos.d1 <- sample(pos.mks[!(pos.mks %in% pos.non)],type.quant[[2]])
  pos.d1.1 <- 1:round(length(pos.d1)/2,0)
  pos.d1.2 <- (round(length(pos.d1)/2,0) + 1):length(pos.d1)
  pos.b3 <- sample(pos.mks[!(pos.mks %in% c(pos.non, pos.d1))],type.quant[[3]])
  pos.d2 <- sample(pos.mks[!(pos.mks %in% c(pos.non, pos.d1, pos.b3))],type.quant[[4]])
  pos.d2.1 <- 1:round(length(pos.d2)/2,0)
  pos.d2.2 <- (round(length(pos.d2)/2,0) + 1):length(pos.d2)
  
  founderfiles <- matrix(NA, nrow = length(mk.names), ncol = 4)
  
  # non-informative
  founderfiles[pos.non, 1:4] <- c(ref.alleles[pos.non], ref.alleles[pos.non], alt.alleles[pos.non], alt.alleles[pos.non])
  # D1.10
  founderfiles[pos.d1[pos.d1.1], 1:4] <- c(ref.alleles[pos.d1[pos.d1.1]], alt.alleles[pos.d1[pos.d1.1]], 
                                           alt.alleles[pos.d1[pos.d1.1]], alt.alleles[pos.d1[pos.d1.1]])
  founderfiles[pos.d1[pos.d1.2], 1:4] <- c(ref.alleles[pos.d1[pos.d1.2]], alt.alleles[pos.d1[pos.d1.2]], 
                                           ref.alleles[pos.d1[pos.d1.2]], ref.alleles[pos.d1[pos.d1.2]])
  # B3.7
  founderfiles[pos.b3, 1:4] <- c(ref.alleles[pos.b3], alt.alleles[pos.b3], ref.alleles[pos.b3], alt.alleles[pos.b3])
  # D2.15
  founderfiles[pos.d2[pos.d2.1], 1:4] <- c(ref.alleles[pos.d2[pos.d2.1]], alt.alleles[pos.d2[pos.d2.1]], 
                                           alt.alleles[pos.d2[pos.d2.1]], alt.alleles[pos.d2[pos.d2.1]])
  founderfiles[pos.d2[pos.d2.2], 1:4] <- c(ref.alleles[pos.d2[pos.d2.2]], alt.alleles[pos.d2[pos.d2.2]], 
                                           ref.alleles[pos.d2[pos.d2.2]], ref.alleles[pos.d2[pos.d2.2]])
  
  # Sample inside genotypes
  founderfiles <- t(apply(founderfiles, 1, function(x) {
    y <- sample(x[1:2])
    z <- sample(x[3:4])
    c(y,z)
  }))
  
  # Sample between parents
  founderfiles <- t(apply(founderfiles, 1, function(x) {
    y <- sample(1:2,1)
    if(y == 2) x[c(3,4,1,2)] else x[-5]
  }
  ))
  
  founderfiles <- cbind(mk.names, founderfiles)
  colnames(founderfiles) <- c("marker", "P1_1", "P1_2", "P2_1", "P2_2")
  
  write.table(founderfiles, file = filename, quote=FALSE, col.names = TRUE, row.names = FALSE, sep = "\t" )
  return(founderfiles)
}

#' Remove inverted SNPs 
#' Sometimes the input linkage map has inversions compared with the physical position, 
#' Here these markers are eliminated
#' @param ref.map data.frame with colunms bp and cM
#' @param thr measure in centimorgan defining the maximum difference accepted for inversions 
#' 
remove_outlier <- function(ref.map, thr=30){
  if(length(which(is.na(ref.map$bp))) > 0)
    ref.map <- ref.map[-which(is.na(ref.map$bp)),]
  
  ref.map$bp <- as.numeric(ref.map$bp)
  
  count <- vector()
  ref.map <- ref.map[order(ref.map$bp),]
  ref.map <- cbind(1:length(ref.map$cM), ref.map)
  while(length(which(diff(ref.map$cM) < -thr)) > 0){
    idx <- which(diff(ref.map$cM) < -thr)
    if(any(idx == 1)) {
      idx <- idx[-which(idx == 1)]
      idx.t <- 1 
    } else { 
      idx.t <- vector()
    }
    test <- ref.map$cM[idx + 1] - ref.map$cM[idx -1]  < -10
    idx.b <- idx[which(test)] + 1
    idx.f <- idx[which(!test)]
    idx.t <- sort(c(idx.t, idx.f, idx.b))
    
    ref.map <- ref.map[-idx.t,]
    count <- c(count , idx.t)
  }
  
  cat("Markers filtered:", length(count), "\n")
  ref.map <- ref.map[,-1]
  return(ref.map)
}


#' Creates map file for PedigreeSim using a reference linkage map or a cM by Mb rate.
#' @param  vcfR.obj object of class vcfR
#' @param ref.map data.frame with two columns nameed cM and bp with numerical values
#' @param cMbyMb numeric value for the cM by Mb rate
create_mapfile <- function(vcfR.obj, ref.map = NULL, cMbyMb = NULL, filename="mapfile.txt", rm_outlier=T, thr=0){
  
  ref.map$cM <- as.numeric(ref.map$cM)
  ref.map$bp <- as.numeric(ref.map$bp)
  pos.vcf <- as.numeric(as.character(vcfR.obj@fix[,2]))
  
  chr = vcfR.obj@fix[,1] 
  pos = vcfR.obj@fix[,2] 
  ref = vcfR.obj@fix[,4] 
  alt = vcfR.obj@fix[,5] 
  
  if(!is.null(cMbyMb)){
    position <- pos.vcf/1000000*cMbyMb
  } else {
    
    # remove inverted
    if(rm_outlier)
      ref.map <- remove_outlier(ref.map, thr)
    
    # remove start and end of chromosome absent in the genetic map
    rm.mks <- unique(c(which(pos.vcf < min(ref.map$bp)), which(pos.vcf > max(ref.map$bp))))
    if(length(rm.mks) > 0){
      pos.vcf <- pos.vcf[-rm.mks]
      chr = chr[-rm.mks] 
      pos = pos[-rm.mks] 
      ref = ref[-rm.mks] 
      alt = alt[-rm.mks] 
    }
    
    # model <- splinefun(ref.map$bp, ref.map$cM, method = "hyman")
    # position <- model(pos.vcf)
    
    # Returns negative values for the gaps
    model <- smooth.spline(ref.map$bp, ref.map$cM)
    position <- predict(model, pos.vcf)$y
    # 
    # model <- gam(ref.map$cM ~ s(ref.map$bp, bs="ps"))
    # position <- predict(model, pos.vcf)
  }
  
  mapfile <- data.frame(marker= paste0(unique(chr), "_", pos.vcf), 
                        chromosome = unique(chr), 
                        position)
  
  write.table(mapfile, file = filename, quote=FALSE, col.names = T, row.names = FALSE, sep = "\t" )
  
  ref_alt_alleles <- data.frame(chr, 
                                pos, 
                                ref, 
                                alt, 
                                pos.map = position,
                                stringsAsFactors = F)
  
  return(list(mapfile, ref_alt_alleles))
}


#' Creates parfile PedigreeSim input
#' By now only available for F1 diploid population
create_parfile <- function(seed, popsize, filename = "parameters.txt"){
  
  parameter <- paste0("PLOIDY = 2
                           MAPFUNCTION = HALDANE
                           MISSING = NA
                           CHROMFILE = chromosome.txt
                           POPTYPE = F1
                           SEED = ", seed,"
                           POPSIZE = ", popsize,"
                           MAPFILE = mapfile.txt
                           FOUNDERFILE = founderfile.txt
                           OUTPUT = sim")
  write.table(parameter, file = filename, quote=FALSE, col.names = FALSE, row.names = FALSE, sep = "\t" )
  return(parameter)
}

#' Creates chrom file PedigreeSim input
#' Only for diploid species
#' @param mapfile output object from create_mapfile
create_chromfile <- function(mapfile, filename = "chromosome.txt"){
  max.pos <- max(mapfile$position)
  chrom <- data.frame("chromosome"= "Chr10", "length"= max.pos, "centromere"=max.pos/2, "prefPairing"= 0.0, "quadrivalents"=0.0)
  write.table(chrom, file= filename, quote = F, col.names = T, row.names = F, sep= "\t")
  return(chrom)
}

#' Codifying phases according with Gusmap
#' for comparisions 
compare_phases <- function(founderfile, ref_alt_alleles, filename = "simulated_phases.txt"){
  founder <- founderfile[,-1]
  simulated_phases <- rep(NA, dim(founder)[1])
  simulated_phases[which(founder[,1] == founder[,3] & founder[,2] == founder[,4])] <- 17 # 1 and 4
  simulated_phases[which(founder[,1] == founder[,4] & founder[,2] == founder[,3])] <- 18 # 2 and 3
  simulated_phases[which(founder[,1] == founder[,3] & founder[,1] == founder[,4] & founder[,1] != founder[,2])] <- 19 # 5 and 8
  simulated_phases[which(founder[,2] == founder[,3] & founder[,2] == founder[,4] & founder[,1] != founder[,2])] <- 20 # 6 and 7
  simulated_phases[which(founder[,1] == founder[,2] & founder[,1] == founder[,3] & founder[,1] != founder[,4])] <- 21 # 9 and 12
  simulated_phases[which(founder[,1] == founder[,2] & founder[,1] == founder[,4] & founder[,1] != founder[,3])] <- 22 # 10 and 11
  simulated_phases[which(founder[,1] == founder[,2] & founder[,1] == founder[,3] & founder[,1] == founder[,4])] <- 23 # 13 and 16
  simulated_phases[which(founder[,1] == founder[,2] & founder[,3] == founder[,4] & founder[,1] != founder[,3])] <- 24 # 14 and 15
  
  simulated_phases <- data.frame(pos=ref_alt_alleles[,2], simulated_phases)
  
  write.table(simulated_phases, file = filename)
  return(simulated_phases)
}
