
# Path
path <- "~/github/errors_workflow/cromwell-executions/F2/0053b719-69f8-466c-b3f0-861f5a7c35d1/"
setwd(path)
cMByMb <- 4.63


# Packages

library(ggplot2)
library(reshape2)
library(vcfR)


method <- c("gatk", "freebayes")
meth.geno <- c("GQ", "polyrad", "updog", "supermassa")

mapfile <- read.table("call-pedsim_files/execution/mapfile.map", header=T)

tot_mks <- read.table("call-pedsim_files/execution/tot_mks.txt")
id <- c(0,1)
p  <- list()

df.diff.dis <- as.list(rep(NA, length(method)))
for(i in 1:length(df.diff.dis)){
  df.diff.dis[[i]] <- as.list(rep(NA, length(meth.geno)))
}

df.per.geno <- df.coverage <- df.diff.tot <- df.diff.median <- df.diff.mean <- vector()

for(i in 1:length(method)){
  ## vcfs
  vcf <- read.vcfR(paste0("call-vcftools_filter/execution/",method[i], ".recode.vcf"))
  gt <- vcf@gt[,-1]
  gt <- melt(gt)
  gt$value <- sapply(strsplit(as.character(gt$value), ":"), "[", 1)
  
  ## Depths
  alt_depth <- read.table(paste0("call-aval_vcf/execution/", method[i], "_alt_depth.txt"), header = T, stringsAsFactors = F)
  alt_depth <- as.data.frame(apply(alt_depth, 2, as.integer))
  alt <- melt(alt_depth)
  
  ref_depth <- read.table(paste0("call-aval_vcf/execution/", method[i], "_ref_depth.txt"), header = T, stringsAsFactors = F)
  ref_depth <- as.data.frame(apply(ref_depth, 2, as.integer))
  ref <- melt(ref_depth)
  df.depths <- data.frame(ind = gt$Var2, gt = gt$value, alt = alt$value, ref=ref$value)
  
  p[[i]] <- ggplot(df.depths, aes(x=ref, y=alt, color=as.factor(gt))) + geom_point(size=0.5) +
    labs(title= paste0(i," depths example"),x="ref", y = "alt", color="Genotypes") 
  
  for(j in 1:length(meth.geno)){
    # Map distance
    maps <- read.table(paste0("call-all_maps/shard-", id[i],"/execution/", method[i],"_map_", meth.geno[j],".txt"), header =T)
    
    poscM <- (as.numeric(as.character(maps[,2]))/1000000)*cMByMb
    poscM.norm <- poscM-poscM[1]
    maps <- cbind(maps, poscM, poscM.norm)
    
    # Difference between distances real and estimated
    diff.dis <- sqrt((maps$poscM.norm - maps$rf)^2)
    df.diff.dis[[i]][[j]] <- assign(paste0(method[i],"_",meth.geno[j]),diff.dis)
    
    # Mean of differences
    diff.dis.mean <- mean(sqrt((maps$poscM.norm - maps$rf)^2))
    df.diff.mean <- c(df.diff.mean, assign(paste0(method[i],"_",meth.geno[j]),diff.dis.mean))
    
    # Median of differences
    diff.dis.median <- median(sqrt((maps$poscM.norm - maps$rf)^2))
    df.diff.median <- c(df.diff.median, assign(paste0(method[i],"_",meth.geno[j]),diff.dis.median))
    
    # Difference of total size
    diff.tot.dis <- sqrt((maps$poscM.norm[length(maps$poscM.norm)] - maps$rf[length(maps$rf)])^2)
    df.diff.tot <- c(df.diff.tot, assign(paste0(method[i],"_",meth.geno[j]),diff.tot.dis))
    
    mapfile.pos <- cbind(mapfile, tot_mks[,2])
    
    # Percentage of the chromosome covered 
    coverage <- maps_gatk$pos[length(maps_gatk$pos)]*100/mapfile.pos[,4][length(mapfile.pos[,4])]
    df.coverage <- c(df.coverage, assign(paste0(method[i],"_",meth.geno[j]),coverage))
    
    # Percentage of rigth genotyping
    per.geno <- read.table(paste0("call-all_maps/shard-",id[i],"/execution/",method[i],"_error_info_",meth.geno[j],".txt"), header =T)
    per.geno <- (sum(per.geno$gabGT == per.geno$methGT)/length(per.geno$methGT)*100)
    df.per.geno <- cbind(df.per.geno, assign(paste0(method[i],"_",meth.geno[j]),per.geno))
    
  }              
}

# Depths
p

# Difference between distances real and estimated
df.diff.dis
# Difference of total size
df.diff.tot
# Median of differences
df.diff.median
# Mean of differences
df.diff.mean

# Percentage of rigth genotyping
df.per.geno
# Percentage of the chromosome covered 
df.coverage
