#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(vcfR)

vcfR.obj <- read.vcfR(args[1])

# Split VCF by potential multiallelic markers combining phased SNPs
split_mult <- function(vcfR.object = NULL, 
                       cross = "outcross", 
                       parent1 = NULL, 
                       parent2 = NULL, 
                       out_file = "position_multi2.txt"){
  
  if (is.null(vcfR.object)) {
    stop("You must specify one vcfR object.")
  }
  if (is.null(parent1) || is.null(parent2)) {
    stop("You must specify samples as parents 1 and 2.")
  }
  if(!is(vcfR.object,"vcfR")){
    stop("You must specify one vcfR object.")
  }
  
  vcf <- vcfR.object
  n.mk <- dim(vcf@gt)[1]
  n.ind <- dim(vcf@gt)[2]-1
  INDS <- dimnames(vcf@gt)[[2]][-1]
  
  MKS <- vcf@fix[,3]
  if (any(MKS == "." | is.na(MKS))) MKS <- paste0(vcf@fix[,1],"_", vcf@fix[,2])
  
  # Geno matrix
  GT_matrix <- matrix(rep(NA,n.ind*n.mk), ncol=n.ind, nrow=n.mk)
  GT <- which(strsplit(vcf@gt[1,1], split=":")[[1]]=="GT")
  
  for(i in 2:(n.ind+1))
    GT_matrix[,i-1] <- unlist(lapply(strsplit(vcf@gt[,i], split=":"), "[[", GT))
  
  CHROM <- vcf@fix[,1]
  POS <- as.numeric(vcf@fix[,2])
  
  # Checking marker segregation according with parents
  P1 <- which(dimnames(vcf@gt)[[2]]==parent1) - 1
  P2 <- which(dimnames(vcf@gt)[[2]]==parent2) - 1
  
  if(length(P1)==0 | length(P2)==0) stop("One or both parents names could not be found in your data")
  
  # This part convert phased genotypes in mnps markers
  if(length(grep("[|]", GT_matrix[,c(P1,P2)])) > 0){
    all_data <- GT_matrix
    all_pos <- POS
    temp_pos <- vector()
    contigs <- unique(CHROM)
    positions <- data.frame()
    # garantee that is the same contig
    for(w in 1:length(contigs)){
      cat(w, "\n")
      idx <- which(CHROM == contigs[w]) 
      GT_matrix <- all_data[idx,]
      POS <- all_pos[idx]
      
      if(length(POS) > 1){
        phased <- grep("[|]", GT_matrix[,P1])
        idx <- which(phased[-1] - phased[-length(phased)] ==1)
        idx.tot <- unique(sort(c(idx, idx +1)))
        idx.p1 <- phased[idx.tot]
        phased <- grep("[|]", GT_matrix[,P2])
        idx <- which(phased[-1] - phased[-length(phased)] ==1)
        idx.tot <- unique(sort(c(idx, idx +1)))
        idx.p2 <- phased[idx.tot]
        idx.tot <- unique(sort(c(idx.p1, idx.p2)))
        
        # Filt NAs unphased heterozygotes
        idx.tot <- idx.tot[-which(grepl(GT_matrix[idx.tot,P1], pattern = "[.]") |  grepl(GT_matrix[idx.tot,P2],pattern = "[.]"))]
        idx.tot <- idx.tot[-which(GT_matrix[idx.tot,P1] == "0/1" |  GT_matrix[idx.tot,P2] == "0/1")]
        idx.tot <- idx.tot[-which(GT_matrix[idx.tot,P1] == "0|0" &  GT_matrix[idx.tot,P2] == "0|0")]
        idx.tot <- idx.tot[-which(GT_matrix[idx.tot,P1] == "1|1" |  GT_matrix[idx.tot,P2] == "1|1")]
        
        idx <- which(idx.tot[-1] - idx.tot[-length(idx.tot)] ==1)
        idx.tot2 <- unique(sort(c(idx, idx +1)))
        idx.tot <- idx.tot[idx.tot2]
        
        #list with haplo
        mnps.num <- split(idx.tot, cumsum(c(1, diff(idx.tot) != 1)))
        pos.mnps <- lapply(mnps.num, function(x) POS[x])
        if(length(unlist(pos.mnps)) > 0){
          positions_temp <- data.frame(chromosome = contigs[[w]], position = unlist(pos.mnps))
        } else {
          positions_temp <- data.frame()
        }
      } else positions_temp <- data.frame()
      positions <- rbind(positions, positions_temp)
    }
    write.table(positions, file = out_file, quote = F, sep = "\t", col.names = F, row.names = F)
  } else {
    positions_null <- data.frame(chromosome=0, position=0)
    write.table(positions_null, file = out_file, quote = F, sep = "\t", col.names = F, row.names = F)
  }
}

split_mult(vcfR.object = vcfR.obj, 
           cross = "outcross", 
           parent1 = args[2], 
           parent2 = args[3], 
           out_file = args[4])
