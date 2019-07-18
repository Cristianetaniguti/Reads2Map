# Packages

library(vcfR)
library(polyRAD)
library(updog)
library(supermassa4onemap)
library(doParallel)
library(onemap)

# Functions

simulations.f2 <- function(i, mean.depth=50, p.mean.depth=50){
  onemap::run_pedsim(chromosome = c("Chr10"),
                     n.marker = c(100), tot.size.cm = c(100), centromere = c(50),
                     n.ind = 200, mk.types = c("A.H.B"), n.types = c(100),
                     pop = "F2", path.pedsim = "~/Programs/PedigreeSim/",
                     name.mapfile = paste0("mapfile",i,".map"),
                     name.founderfile=paste0("founderfile.",i,".gen"),
                     name.chromfile=paste0("sim",i,".chrom"), name.parfile=paste0("sim",i,".par"),
                     name.out=paste0("sim_out.f2.",i))
  
  onemap::pedsim2vcf(inputfile=paste0("sim_out.f2.",i,"_genotypes.dat"),
                     map.file=paste0("mapfile",i,".map"),
                     chrom.file=paste0("sim",i,".chrom"),
                     out.file=paste0("out.f2.",i,".vcf"),
                     miss.perc = 0,
                     counts=TRUE,
                     mean.depth= mean.depth,
                     p.mean.depth = p.mean.depth,
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
}

map.build <- function(onemap.obj, cores=6){
  cl<-makeCluster(cores)
  twopts <- parLapply(cl, onemap.obj,
                      function(x) onemap::rf_2pts(x))
  
  for(w in 1:length(twopts)){
    twopts[[w]]$data.name <- paste0("data.name.",w)
    assign(paste0("twopts.", w), twopts[[w]])
  }
  
  for (i in 1:n.sim){assign(paste0("twopts.",i), twopts[[i]])}
  twopts.names <- paste0("twopts.",1:n.sim)
  
  seq1 <- parLapply(cl, twopts,
                    function(x) onemap::make_seq(x, "all"))
  
  for(i in 1:n.sim){
    seq1[[i]]$data.name <- paste0("data.name.",i)
    seq1[[i]]$twopt <- paste0("twopts.",i)
  }
  
  clusterExport(cl, c(data.names, twopts.names), envir=environment())
  
  maps <- parLapply(cl, seq1,
                    function(x) onemap::map(x))
  
  stopCluster(cl)
  dist <- dist.ind <- list()
  for(i in 1:n.sim){
    dist.ind[[i]] <- round(haldane(maps[[i]]$seq.rf),3)
    dist[[i]] <- c(0,round(cumsum(haldane(maps[[i]]$seq.rf)),3))
  }
  
  return(list(maps, dist.ind, dist))
}



# F2 simulations

depth <- c(50,10,2)
for(j in 1:length(depth)){
  
  cores <- 6
  n.sim <- 100
  idx <- as.list(1:n.sim)
  
  cl<-makeCluster(cores)
  clusterExport(cl, c("simulations.f2", "depth", "j"), envir=environment())
  parLapply(cl, idx, function(x) simulations.f2(i=x, mean.depth = depth[j], p.mean.depth = depth[j]))
  stopCluster(cl)
  
  
  files <- list()
  for(i in 1:n.sim)
    files[[i]] <- paste0("out.f2.",i, ".vcf")
  
  cl<-makeCluster(cores)
  
  data.vcf <- parLapply(cl, files,
                        function(x) vcfR::read.vcfR(file=x))
  onemap.vcf <- parLapply(cl,
                          data.vcf,
                          function(x) onemap::onemap_read_vcfR(vcfR.object = x, 
                                                               cross = "f2 intercross", 
                                                               parent1 = "P1", 
                                                               parent2 = "P2", 
                                                               f1="F1"))
  stopCluster(cl)
  
  ## Default
  
  for (i in 1:n.sim){assign(paste("data.name.",i, sep=""), onemap.vcf[[i]])}
  data.names <- paste0("data.name.",1:n.sim)
  
  map.df <- map.build(onemap.vcf)
  
  
  ## GQ
  # Here there are not the GQ information but at the workflow this error will also be considered
  
  ## updog
  
  cl<-makeCluster(cores)
  clusterExport(cl, c("data.vcf", "onemap.vcf", "flexdog"), envir=environment())
  
  onemap.updog <- parLapply(cl, idx,
                            function(x)
                              onemap::updog_error(
                                vcfR.object=data.vcf[[x]],
                                onemap.object = onemap.vcf[[x]],
                                vcf.par = "AD",
                                parent1 = "P1",
                                parent2 = "P2",
                                f1="F1",
                                recovering = TRUE,
                                mean_phred = 20,
                                cores = 6,
                                depths = NULL))
  stopCluster(cl)
  
  for (i in 1:n.sim){assign(paste("data.name.",i, sep=""), onemap.updog[[i]])}
  data.names <- paste0("data.name.",1:n.sim)
  
  map.up <- map.build(onemap.updog)
  
  ## polyRAD
  
  cl<-makeCluster(cores)
  clusterExport(cl, c("files", "onemap.vcf"), envir=environment())
  poly.onemap <- parLapply(cl, idx,
                           function(x)
                             onemap::polyRAD_error(vcf=files[[x]], 
                                                   onemap.obj = onemap.vcf[[x]],
                                                   parent1="P1",
                                                   parent2="P2",
                                                   f1="F1",
                                                   crosstype="f2 intercross"))
  
  
  stopCluster(cl)
  
  map.poly <- map.build(poly.onemap)
  
  ## supermassa
  cl<-makeCluster(cores)
  clusterExport(cl, c("data.vcf", "onemap.vcf"), envir=environment())
  super.onemap <- parLapply(cl, idx,
                            function(x)
                              supermassa4onemap::supermassa_error(vcfR.object=data.vcf[[x]],
                                                                  onemap.object = onemap.vcf[[x]],
                                                                  vcf.par = "AD",
                                                                  parent1 = "P1",
                                                                  parent2 = "P2",
                                                                  f1="F1",
                                                                  recovering = TRUE,
                                                                  mean_phred = 20,
                                                                  cores = 3,
                                                                  depths = NULL))
  
  stopCluster(cl)
  
  map.super <- map.build(super.onemap)
  
  # Graphics
  
  maps.names <- c("map.df",
                  "map.up",
                  "map.poly",
                  "map.super")
  
  # Total size
  
  genotype <- rep(maps.names, each=n.sim)
  
  tot.size <- list()
  for(i in 1:length(maps.names))
    tot.size[[i]] <- sapply(get(maps.names[i])[[3]], function(x) x[length(x)])
  tot.size <- unlist(tot.size)
  
  tot.size.df <- data.frame(genotype, tot.size)
  
  p <- ggplot(tot.size.df, aes(x=genotype, y=tot.size)) + geom_boxplot() +
    labs(title= "Total size" ,x="genotype and probs", y = "cM")   + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          panel.background = element_rect(fill = "lightblue",
                                          colour = "lightblue",
                                          size = 0.5, linetype = "solid")) +
    geom_hline(yintercept = 12, color = "darkgreen") 
  
  
  
  ggsave(p,filename = paste0("Tot.size.jpg"))
  p
  
  # Distances between markers
  
  df.dist.part <- data.frame()
  for(i in 1:length(maps.names)){
    dist.part <- unlist(get(maps.names[i])[[2]])
    df.dist.part <- rbind(df.dist.part, data.frame(geno=rep(maps.names[i], length(dist.part)), dist = dist.part))
  }
  
  
  p <- ggplot(df.dist.part, aes(x=geno, y=dist)) + geom_boxplot() +
    labs(title= "Distances" ,x="genotype and probs", y = "cM")   + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          panel.background = element_rect(fill = "lightblue",
                                          colour = "lightblue",
                                          size = 0.5, linetype = "solid")) +
    geom_hline(yintercept = 1, color = "darkgreen") 
  
  ggsave(p,filename = paste0("Dists.jpg"))
  p
}

##################
# Outcrossing
##################

cores <- 2
n.sim <- 2
idx <- as.list(1:n.sim)

simulations.out <- function(i){
  onemap::run_pedsim(chromosome = c("Chr10"),
                     n.marker = c(100), tot.size.cm = c(100), centromere = c(50),
                     n.ind = 200, mk.types = c("B3.7", "D1.10", "D2.15"), n.types = c(30,35,35),
                     pop = "F1", path.pedsim = "~/Programs/PedigreeSim/",
                     name.mapfile = "mapfile.map",
                     name.founderfile=paste0("founderfile.",i,".gen"),
                     name.chromfile="sim.chrom", name.parfile="sim.par",
                     name.out=paste0("sim_out.",i))
  
  onemap::pedsim2vcf(inputfile=paste0("sim_out.f2.",i,"_genotypes.dat"),
                     map.file="mapfile.map",
                     chrom.file="sim.chrom",
                     out.file=paste0("out.",i,".vcf"),
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
}

cl<-makeCluster(cores)
clusterExport(cl, c("simulations.out"), envir=environment())
parLapply(cl, idx, function(x) simulations.out(i=x))
stopCluster(cl)

# use founderfiles to save the phases in outcrossing

real.phase <- read.table("founderfile.1.gen", header=T)
opgp.real <- rep(NA, dim(real.phase)[1])

index <- which(real.phase[,2] == "a" & real.phase[,3] == "b")

opgp.real[index][which(real.phase[index, 4] == "a" & real.phase[index, 5] == "b")] <- 1

opgp.real[index][which(real.phase[index, 4] == "b" & real.phase[index, 5] == "a")] <- 3

opgp.real[index][which(real.phase[index, 4] == "a" & real.phase[index, 5] == "a")] <- 5

opgp.real[index][which(real.phase[index, 4] == "b" & real.phase[index, 5] == "b")] <- 7


index <- which(real.phase[,2] == "b" & real.phase[,3] == "a")

opgp.real[index][which(real.phase[index, 4] == "a" & real.phase[index, 5] == "b")] <- 2

opgp.real[index][which(real.phase[index, 4] == "b" & real.phase[index, 5] == "a")] <- 4

opgp.real[index][which(real.phase[index, 4] == "a" & real.phase[index, 5] == "a")] <- 6

opgp.real[index][which(real.phase[index, 4] == "b" & real.phase[index, 5] == "b")] <- 8


index <- which(real.phase[,2] == "a" & real.phase[,3] == "a")

opgp.real[index][which(real.phase[index, 4] == "a" & real.phase[index, 5] == "b")] <- 9

opgp.real[index][which(real.phase[index, 4] == "b" & real.phase[index, 5] == "a")] <- 10

opgp.real[index][which(real.phase[index, 4] == "a" & real.phase[index, 5] == "a")] <- 13

opgp.real[index][which(real.phase[index, 4] == "b" & real.phase[index, 5] == "b")] <- 14


index <- which(real.phase[,2] == "b" & real.phase[,3] == "b")

opgp.real[index][which(real.phase[index, 4] == "a" & real.phase[index, 5] == "b")] <- 11

opgp.real[index][which(real.phase[index, 4] == "b" & real.phase[index, 5] == "a")] <- 12

opgp.real[index][which(real.phase[index, 4] == "a" & real.phase[index, 5] == "a")] <- 15

opgp.real[index][which(real.phase[index, 4] == "b" & real.phase[index, 5] == "b")] <- 16

# For onemap phases 1 and 4; 2 and 3 are the same, so

opgp.real[which(opgp.real == 1 | opgp.real == 4)] <- 17
opgp.real[which(opgp.real == 2 | opgp.real == 3)] <- 18

phaseToOPGP  <- function(x){
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
  if (class(get(x$data.name, pos = 1))[2] == "outcross") {
    link.phases <- apply(link.phases, 1, function(x) paste(as.character(x), collapse = "."))
    parents <- matrix("", length(x$seq.num), 4)
    for (i in 1:length(x$seq.num))
      parents[i, ] <- onemap::return <- geno(get(x$data.name, pos = 1)$segr.type[x$seq.num[i]], link.phases[i])
    ## Our code below
    #transpose the the parents and set to baseline
    parents[which(parents == 'a')] <-'A'
    parents[which(parents == 'b')] <- 'B'
    
    parents = t(parents)
    if(parents[1,which(apply(parents[1:2,],2,function(x) !(all(x=='A'))))[1]] == 'B')
      parents[1:2,] <- parents[2:1,]
    if(parents[3,which(apply(parents[3:4,],2,function(x) !(all(x=='A'))))[1]] == 'B')
      parents[3:4,] <- parents[4:3,]
    
    ## Now from the parental haplotypes, determine the OPGPs
    return(GUSMap:::parHapToOPGP(parents))
  }
}



map.build <- function(onemap.obj, cores=4){
  cl<-makeCluster(cores)
  twopts <- parLapply(cl, onemap.obj,
                      function(x) onemap::rf_2pts(x))
  
  for(w in 1:length(twopts)){
    twopts[[w]]$data.name <- paste0("data.name.",w)
    assign(paste0("twopts.", w), twopts[[w]])
  }
  
  for (i in 1:n.sim){assign(paste0("twopts.",i), twopts[[i]])}
  twopts.names <- paste0("twopts.",1:n.sim)
  
  seq1 <- parLapply(cl, twopts,
                    function(x) onemap::make <- seq(x, "all"))
  
  for(i in 1:n.sim){
    seq1[[i]]$data.name <- paste0("data.name.",i)
    seq1[[i]]$twopt <- paste0("twopts.",i)
  }
  
  clusterExport(cl, c(data.names, twopts.names, "phaseToOPGP_OM", "return_geno"), envir=environment())
  
  maps <- parLapply(cl, seq1,
                    function(x) onemap::map(x))
  
  phases <- parLapply(cl, maps,
                      function(x) phaseToOPGP <- OM(x))
  
  stopCluster(cl)
  perc.phase <- dist <- dist.ind <- list()
  for(i in 1:n.sim){
    names(phases[[i]]) <- colnames(get(paste0("data.name.",i))$geno)
    phases[[i]][which(phases[[i]] == 1 | phases[[i]] == 4)] <- 17
    phases[[i]][which(phases[[i]] == 2 | phases[[i]] == 3)] <- 18
    
    pos1 <- match(names(phases[[i]]), names(real.phase))
    perc.phase[[i]] <- 100*sum(phases[[i]] == real.phase[pos1])/length(real.phase)
    
    dist.ind[[i]] <- round(haldane(maps[[i]]$seq.rf),3)
    dist[[i]] <- c(0,round(cumsum(haldane(maps[[i]]$seq.rf)),3))
    names(dist[[i]]) <- names(phases[[i]])
  }
  
  return(list(maps, dist.ind, dist, phases, perc.phase))
}

files <- list()
for(i in 1:n.sim)
  files[[i]] <- paste0("out.f2.",i, ".vcf")

cl<-makeCluster(cores)

data.vcf <- parLapply(cl, files,
                      function(x) vcfR::read.vcfR(file=x))

onemap.vcf <- parLapply(cl,
                        data.vcf,
                        function(x) onemap::onemap_read_vcfR(vcfR.object = x, cross = "outcross", parent1 = "P1", parent2 = "P2"))
stopCluster(cl)

### ta dando pau aqui, na leitura de outcross

# Here there are not the GQ information

# Default

for (i in 1:n.sim){assign(paste("data.name.",i, sep=""), onemap.vcf[[i]])}
data.names <- paste0("data.name.",1:n.sim)

map.df <- map.build(onemap.vcf)


## updog

cl<-makeCluster(cores)
clusterExport(cl, c("data.vcf", "onemap.vcf", "flexdog"), envir=environment())

onemap.updog <- parLapply(cl, idx,
                          function(x)
                            onemap::updog_error(
                              vcfR.object=data.vcf[[x]],
                              onemap.object = onemap.vcf[[x]],
                              vcf.par = "AD",
                              parent1 = "P1",
                              parent2 = "P2",
                              f1="F1",
                              recovering = TRUE,
                              mean_phred = 20,
                              cores = 6,
                              depths = NULL))
stopCluster(cl)

for (i in 1:n.sim){assign(paste("data.name.",i, sep=""), onemap.updog[[i]])}
data.names <- paste0("data.name.",1:n.sim)

map.up <- map.build(onemap.updog)

