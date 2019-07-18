###########
# Packages
###########

library(vcfR)
library(polyRAD)
library(updog)
library(supermassa4onemap)
library(doParallel)

####################################
# Testing for outcrossing species
####################################

###########
# Using GQ
###########

vcf.out <- read.vcfR("vcf_example_out.vcf")
out <- onemap_read_vcfR(vcfR.object = vcf.out, cross = "outcross", 
                        parent1 = "P1", parent2 = "P2")

segr <- test_segregation(out)
plot(segr)

twopts <- rf_2pts(out)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg1  <- make_seq(lgs,2)
map.lg1.1 <- map(lg1)

genotypes_errors <- extract_depth(vcfR.object = vcf.out, onemap.object = out, vcf.par = "GQ", parent1 = "P1", 
                                 parent2 = "P2", mean_phred = 20, recovering = FALSE)

new.errors <- create_probs(out, 
                           genotypes_errors = genotypes_errors, 
                           global_error = NULL, 
                           genotypes_probs = NULL)

twopts <- rf_2pts(new.errors)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg1  <- make_seq(lgs,2)
map.lg1.1 <- map(lg1)

########
# updog
########

old.updog <- updog_error(vcfR.object = vcf.out, onemap.object = out, vcf.par = "AD", parent1 = "P1", 
                         parent2 = "P2", recovering = TRUE, mean_phred = 20, cores = 3, 
                         depths = NULL)

segr <- test_segregation(old.updog)
plot(segr)

twopts <- rf_2pts(old.updog)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg1 <- make_seq(lgs,2)
map.lg1 <- map(lg1)

##########
# PolyRAD
##########

poly.test <- VCF2RADdata("vcf_example_out.vcf", phaseSNPs = FALSE, 
                         min.ind.with.reads = 0,
                         min.ind.with.minor.allele = 0)


poly.test <- SetDonorParent(poly.test, "P1")
poly.test <- SetRecurrentParent(poly.test, "P2")

mydata2 <- PipelineMapping2Parents(poly.test, 
                                   freqAllowedDeviation = 0.06,
                                   useLinkage = FALSE,
                                   minLikelihoodRatio = 2)

Export_MAPpoly(mydata2, "test")

genotypes <- read.table("test", skip=12)

pos <- sapply(strsplit(as.character(genotypes$V1), split = "_"),"[",1)

pos.onemap <- colnames(out$geno)
genotypes <- genotypes[which(pos%in%pos.onemap),]
keep.mks <- which(pos.onemap%in%pos)

# Atualizar geno
out$geno <- out$geno[,keep.mks]

new.geno <- vector()
for(i in 1:dim(genotypes)[1]){
  if(which.max(genotypes[i,3:5]) == 3){
    new.geno[i] <- 3
  }else if(which.max(genotypes[i,3:5]) == 2){
    new.geno[i] <- 2
  } else if(which.max(genotypes[i,3:5]) == 1){
    new.geno[i] <- 1
  }
}

new.geno <- matrix(new.geno,nrow = out$n.ind, ncol = length(keep.mks) )
colnames(new.geno) <- colnames(out$geno)
rownames(new.geno) <- rownames(out$geno)

# Mudando a ordem
genotypes <- genotypes[order(genotypes$V2),]

# Quanto mudou
sum(new.geno == out$geno)/length(new.geno)

out$geno <- new.geno
# Remover marcadores
out$n.mar <- length(keep.mks)
out$segr.type <- out$segr.type[keep.mks]
out$segr.type.num <- out$segr.type.num[keep.mks]
out$CHROM <- out$CHROM[keep.mks]
out$POS <- out$POS[keep.mks]

polyrad.one <- create_probs(onemap.obj = out, genotypes_probs =  genotypes[,3:5])
head(polyrad.one$error)

twopts <- rf_2pts(polyrad.one)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,2)
map.lg2 <- map(lg2)


####################################
# Testing for f2 populations
####################################
# The example was without the F1
fix.f2 <- read.table("vcf_example_f2.vcf")
head(fix.f2)
new.f2 <- cbind(fix.f2[,1:11], rep("0/1:7,4:11:99:111,0,219",length(fix.f2[,1])), fix.f2[,12:dim(fix.f2)[2]])
write.table(new.f2, file = "new.f2", quote = F, sep = "\t", col.names = F, row.names = F)

f2.vcf <- read.vcfR("vcf_example_f2.new.vcf")

###########
# Using GQ
###########

f2 <- onemap_read_vcfR(vcfR.object = f2.vcf, cross = "f2 intercross", parent1 = "P1", parent2 = "P2", f1 = "F1")

segr <- test_segregation(f2)
plot(segr)

twopts <- rf_2pts(f2)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg1  <- make_seq(lgs,2)
map.lg1.1 <- map(lg1)

genotypes_errors <- extract_depth(vcfR.object = f2.vcf, onemap.object = f2, vcf.par = "GQ", parent1 = "P1", 
                                  parent2 = "P2", f1 = "F1", mean_phred = 20, recovering = FALSE)

new.errors <- create_probs(f2, 
                           genotypes_errors = genotypes_errors, 
                           global_error = NULL, 
                           genotypes_probs = NULL)

twopts <- rf_2pts(new.errors)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg1  <- make_seq(lgs,2)
map.lg1.1 <- map(lg1)

########
# updog
########

old.updog <- updog_error(vcfR.object = f2.vcf, onemap.object = f2, vcf.par = "AD", parent1 = "P1", 
                         parent2 = "P2", recovering = TRUE, mean_phred = 20, cores = 3, 
                         depths = NULL)

segr <- test_segregation(old.updog)
plot(segr)

twopts <- rf_2pts(old.updog)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg1 <- make_seq(lgs,2)
map.lg1 <- map(lg1)

##########
# PolyRAD
##########

poly.test <- VCF2RADdata("vcf_example_f2.new.vcf", phaseSNPs = FALSE, 
                         min.ind.with.reads = 0,
                         min.ind.with.minor.allele = 0)


poly.test <- SetDonorParent(poly.test, "F1")
poly.test <- SetRecurrentParent(poly.test, "F1")

mydata2 <- PipelineMapping2Parents(poly.test, 
                                   freqAllowedDeviation = 0.06,
                                   useLinkage = FALSE)

Export_MAPpoly(mydata2, "test")

genotypes <- read.table("test", skip=12)

pos <- sapply(strsplit(as.character(genotypes$V1), split = "_"),"[",1)

pos.onemap <- colnames(f2$geno)
genotypes <- genotypes[which(pos%in%pos.onemap),]
keep.mks <- which(pos.onemap%in%pos)

# Remove parents
genotypes <- genotypes[-which(genotypes[,2]%in%"P1"),]
genotypes <- genotypes[-which(genotypes[,2] %in%"P2"),]

# Atualizar geno
f2$geno <- f2$geno[,keep.mks]

new.geno <- vector()
for(i in 1:dim(genotypes)[1]){
  if(which.max(genotypes[i,3:5]) == 3){
    new.geno[i] <- 3
  }else if(which.max(genotypes[i,3:5]) == 2){
    new.geno[i] <- 2
  } else if(which.max(genotypes[i,3:5]) == 1){
    new.geno[i] <- 1
  }
}

new.geno <- matrix(new.geno,nrow = f2$n.ind, ncol = length(keep.mks))
colnames(new.geno) <- colnames(f2$geno)
rownames(new.geno) <- rownames(f2$geno)

# Mudando a ordem
genotypes <- genotypes[order(genotypes$V2),]

# Quanto mudou
sum(new.geno == f2$geno)/length(new.geno)

f2$geno <- new.geno

# Remover marcadores
f2$n.mar <- length(keep.mks)
f2$segr.type <- f2$segr.type[keep.mks]
f2$segr.type.num <- f2$segr.type.num[keep.mks]
f2$CHROM <- f2$CHROM[keep.mks]
f2$POS <- f2$POS[keep.mks]

polyrad.one <- create_probs(onemap.obj = f2, genotypes_probs =  genotypes[,3:5])
head(polyrad.one$error)

twopts <- rf_2pts(polyrad.one)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,2)
map.lg2 <- map(lg2)

######################################
# Testing with GATK simulated results
######################################

###########
# Using GQ
###########

f2.vcf <- read.vcfR("family1_gatk.vcf")

f2 <- onemap_read_vcfR(vcfR.object = f2.vcf, cross = "f2 intercross", parent1 = "P1", parent2 = "P2", f1 = "F1")

segr <- test_segregation(f2)
plot(segr)

twopts <- rf_2pts(f2)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg1  <- make_seq(lgs,1)
map.lg1.1 <- map(lg1)

genotypes_errors <- extract_depth(vcfR.object = f2.vcf, onemap.object = f2, vcf.par = "GQ", parent1 = "P1", 
                                  parent2 = "P2", f1 = "F1", mean_phred = 20, recovering = FALSE)

new.errors <- create_probs(f2, 
                           genotypes_errors = genotypes_errors, 
                           global_error = NULL, 
                           genotypes_probs = NULL)

twopts <- rf_2pts(new.errors)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg1  <- make_seq(lgs,1)
map.lg1.1 <- map(lg1)

########
# updog
########

old.updog <- updog_error(vcfR.object = f2.vcf, onemap.object = f2, vcf.par = "AD", parent1 = "P1", 
                         parent2 = "P2",f1 = "F1", recovering = TRUE, mean_phred = 20, cores = 3, 
                         depths = NULL)

segr <- test_segregation(old.updog) # Verificar pq o teste ta zuado!
plot(segr)

twopts <- rf_2pts(old.updog)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg1 <- make_seq(lgs,1)
map.lg1 <- map(lg1)

##########
# PolyRAD
##########

poly.onemap <- polyRAD_error(vcf="family1_gatk.vcf", 
                             onemap.obj = f2,
                             parent1="P1",
                             parent2="P2",
                             f1="F1",
                             crosstype="f2 intercross")


twopts <- rf_2pts(poly.onemap)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,2)
map.lg1.poly <- map(lg2)

# fazer teste com exemplos

data("vcf_example_f2")
segr <- test_segregation(vcf_example_f2)
plot(segr)
test.df <- create_probs(vcf_example_f2)
twopts <- rf_2pts(test.df)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,2)
map.lg2 <- map(lg2)

data("vcf_example_riself")
segr <- test_segregation(vcf_example_riself)
plot(segr)

test.df <- create_probs(vcf_example_riself)
twopts <- rf_2pts(test.df)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,2)
map.lg2 <- map(lg2)

data("vcf_example_f2")
test.df <- create_probs(vcf_example_f2)
twopts <- rf_2pts(test.df)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,2)
map.lg2 <- map(lg2)

data("vcf_example_bc")
test.df <- create_probs(vcf_example_bc)
twopts <- rf_2pts(test.df)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,2)
map.lg2 <- map(lg2)


###############################
# Simulations to test the code
###############################

run_pedsim(chromosome = c("Chr10"), n.marker = c(100), tot.size.cm = c(100), centromere = c(50),
           n.ind = 200, mk.types = c("A.H.B"), n.types = c(100), pop = "F2", path.pedsim = "~/Programs/PedigreeSim/",
           name.mapfile = "mapfile.map", name.founderfile="founderfile.gen", name.chromfile="sim.chrom", name.parfile="sim.par",
           name.out="sim_out.f2")

pedsim2vcf(inputfile="sim_out.f2_genotypes.dat", 
           map.file="mapfile.map", 
           chrom.file="sim.chrom",
           out.file="out.f2.vcf", 
           miss.perc = 0, 
           counts=TRUE, 
           mean.depth=50, 
           p.mean.depth = 50, 
           disper.par=2, 
           chr.mb= 10, 
           method = "updog", 
           mean.phred=20, 
           bias=1, 
           od=0.001,
           pos=NULL,
           haplo.ref=NULL,
           chr=NULL,
           phase = FALSE)


f2.vcf <- read.vcfR("out.f2.vcf")

f2 <- onemap_read_vcfR(vcfR.object = f2.vcf, cross = "f2 intercross", parent1 = "P1", parent2 = "P2", f1 = "F1")

########
# updog
########

old.updog <- updog_error(vcfR.object = f2.vcf, onemap.object = f2, vcf.par = "AD", parent1 = "P1", 
                         parent2 = "P2",f1 = "F1", recovering = TRUE, mean_phred = 20, cores = 3, 
                         depths = NULL)

segr <- test_segregation(old.updog) 
plot(segr)

twopts <- rf_2pts(old.updog)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg1 <- make_seq(lgs,1)
map.lg1 <- map(lg1)

##########
# PolyRAD
##########

poly.test <- VCF2RADdata("out.f2.vcf", phaseSNPs = FALSE, 
                         min.ind.with.reads = 0,
                         min.ind.with.minor.allele = 0)


poly.test <- SetDonorParent(poly.test, "F1")
poly.test <- SetRecurrentParent(poly.test, "F1")

mydata2 <- PipelineMapping2Parents(poly.test, 
                                   freqAllowedDeviation = 0.06,
                                   useLinkage = FALSE)

Export_MAPpoly(mydata2, "test")

genotypes <- read.table("test", skip=12)

any(is.na(genotypes[,3:5]))

pos <- sapply(strsplit(as.character(genotypes$V1), split = "_"),"[",1)

pos.onemap <- colnames(f2$geno)
genotypes <- genotypes[which(pos%in%pos.onemap),]
keep.mks <- which(pos.onemap%in%pos)

# Remove parents
genotypes <- genotypes[-which(genotypes[,2]%in%"P1"),]
genotypes <- genotypes[-which(genotypes[,2] %in%"P2"),]

# Atualizar geno
f2$geno <- f2$geno[,keep.mks]

new.geno <- vector()
for(i in 1:dim(genotypes)[1]){
  if(which.max(genotypes[i,3:5]) == 3){
    new.geno[i] <- 3
  }else if(which.max(genotypes[i,3:5]) == 2){
    new.geno[i] <- 2
  } else if(which.max(genotypes[i,3:5]) == 1){
    new.geno[i] <- 1
  }
}

new.geno <- matrix(new.geno,nrow = f2$n.ind, ncol = length(keep.mks))
colnames(new.geno) <- colnames(f2$geno)
rownames(new.geno) <- rownames(f2$geno)

# Mudando a ordem
genotypes <- genotypes[order(genotypes$V2),]

# Quanto mudou
1- sum(new.geno == f2$geno)/length(new.geno)

f2$geno <- new.geno

# Remover marcadores
f2$n.mar <- length(keep.mks)
f2$segr.type <- f2$segr.type[keep.mks]
f2$segr.type.num <- f2$segr.type.num[keep.mks]
f2$CHROM <- f2$CHROM[keep.mks]
f2$POS <- f2$POS[keep.mks]

polyrad.one <- create_probs(onemap.obj = f2, genotypes_probs =  genotypes[,3:5])
head(polyrad.one$error)

twopts <- rf_2pts(polyrad.one)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,1)
map.lg1.poly <- map(lg2)

#############
# Supermassa
#############

f2.vcf <- read.vcfR("out.f2.vcf")

f2 <- onemap_read_vcfR(vcfR.object = f2.vcf, cross = "f2 intercross", parent1 = "P1", parent2 = "P2", f1 = "F1")


supermassa.aval <- supermassa_error(vcfR.object=f2.vcf,
                                    onemap.object = f2,
                                    vcf.par = "AD",
                                    parent1 = "P1",
                                    parent2 = "P2",
                                    f1="F1",
                                    recovering = TRUE,
                                    mean_phred = 20,
                                    cores = 3,
                                    depths = NULL)

supermassa.aval <- create_probs(supermassa.aval, genotypes_errors = supermassa.aval$error)

twopts <- rf_2pts(supermassa.aval)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,1)
map.lg1.super <- map(lg2)

