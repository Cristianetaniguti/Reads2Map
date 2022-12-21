version 1.0

# If there is no reference VCF
task CreatePedigreeSimulatorInputs {
  input {
    File snps
    File indels
    Float? cmBymb
    File ref
    File ref_fai
    Int seed
    Int popsize
    Int ploidy # Only ploidy 2 available
    File? doses
    String cross
  }

  Int disk_size = ceil(size(snps, "GiB") * 2 + size(indels, "GiB") * 2 + size(ref, "GiB") + size(doses, "GiB") * 2)
  Int memory_size = 4000

  command <<<

      R --vanilla --no-save <<RSCRIPT
      library(Reads2MapTools)
      # Needs optimization
      ploidy <- as.integer("~{ploidy}")
      doses <- read.table("~{doses}", sep=",")
      snps <- read.table("~{snps}", stringsAsFactors = FALSE)
      colnames(snps) <- c("chr", "pos.befor", "pos.after", "ref", "alt")
      indels <- read.table("~{indels}", stringsAsFactors = FALSE)
      colnames(indels) <- c("chr", "pos.befor", "pos.after", "sinal", "n.bases", "bases")
      pos.ref <- indels[,2]
      sinal <- indels[,4]


      # Pirs output do not show the last base before the polimorphism
      # When it is a negative indel compared to the reference, the last position before the pointed
      # is the last base before the polimorphism
      pos.ref[which(sinal=="-")] <- pos.ref[which(sinal=="-")] -1

      # search last base before the indels (information needed by VCF)
      int <- paste0(indels[1,1],":",pos.ref,"-",pos.ref)
      sep.idx <- as.integer(length(int)/1000)
      sep <- rep(1:sep.idx, each =1000)
      sep <- c(sep, rep(sep.idx+1, each=length(int)-length(sep)))
      sep <- split(int, sep)

      bases <- list()
      for(i in 1:length(sep))
        bases[[i]] <- system(paste(paste("samtools faidx", "~{ref}"), paste(sep[[i]], collapse =  " ")), intern = T)

      bases <- do.call(c, bases)

      bases.bf <- matrix(bases, ncol=2, byrow = T)[,2]
      alt <- bases.bf
      tmp <- paste0(bases.bf[which(sinal=="+")], indels[,6][which(sinal=="+")])
      alt[which(sinal=="+")] <- tmp
      ref <- bases.bf
      tmp <- paste0(bases.bf[which(sinal=="-")], indels[,6][which(sinal=="-")])
      ref[which(sinal=="-")] <- tmp

      # the position in the vcf and in the map are according with the reference genome
      ref_alt_alleles <- data.frame(chr = c(snps[,1], indels[,1]), pos = c(snps[,2], pos.ref),
                            ref = c(snps[,4], ref), alt = c(snps[,5],alt), stringsAsFactors = F)


      ref_alt_alleles <- ref_alt_alleles[order(ref_alt_alleles[,2]),]
      n.marker <- dim(ref_alt_alleles)[1]

      ## Map file
      # Marker names
      marker1 <- "M"
      marker2 <- 1:n.marker
      marker2 <- sprintf("%03d", marker2)
      marker <-paste0(marker1,marker2)

      # Chromossome and position
      pos.map <- (ref_alt_alleles[,2]/1000000) * ~{cmBymb}
      map_file <- data.frame(marker=marker, chromosome=ref_alt_alleles[,1], position= pos.map)
      write.table(map_file, file = paste0("mapfile.txt"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

      ref_alt_alleles <- cbind(ref_alt_alleles, pos.map)
      write.table(ref_alt_alleles, file="ref_alt_alleles.txt")

      ## Founderfile
      ref.alleles <- ref_alt_alleles[,3]
      alt.alleles  <- ref_alt_alleles[,4]

      if("~{cross}" == "F1"){
        doses <- doses[c(1,length(doses),2:(length(doses)-1))]
        ploidys <- c(0:ploidy)
        ploidys <- ploidys[c(1,length(doses),2:(length(doses)-1))]

        founder1.df <- matrix(NA, ncol = ploidy, nrow = length(ref.alleles))

        founder2.df <- matrix(NA, ncol = ploidy, nrow = length(ref.alleles))

        idx <- 1:length(ref.alleles)
        for(i in 1:length(doses)){
          size <- round((doses[i]/100)*length(ref.alleles))
          if(i == 1){ # homozigote for reference
            idx.both <- sample(idx, as.numeric(size)*2) # It will not have monomorphic markers
            idx.p1 <- idx.both[1:as.numeric(size)]
            idx.p2 <- idx.both[(as.numeric(size)+1):(as.numeric(size)*2)]
            founder1.df[idx.p1,] <- ref.alleles[idx.p1]
            founder2.df[idx.p2,] <- ref.alleles[idx.p2]
            idx.p1.tot <- idx[-idx.p1]
            idx.p2.tot <- idx[-idx.p2]
          } else if(i == 2){ # homozigote for alternative
            idx.p1 <- sample(idx.p1.tot, as.numeric(size)) # select remaining lines for P1
            idx.p2 <- vector() # select remaining lines for P2
            for(w in 1:(as.numeric(size))){
              idx.p2[w] <- sample(idx.p2.tot, 1)
              while(any(idx.p1 %in% idx.p2[w])){ # the line selected in P2 can not be the same of P1
                idx.p2[w] <- sample(idx.p2.tot, 1)
              }
              idx.p2.tot <- idx.p2.tot[-which(idx.p2.tot%in%idx.p2)]
            }
            idx.p1.tot <- idx.p1.tot[-which(idx.p1.tot%in%idx.p1)]
            founder1.df[idx.p1,] <- alt.alleles[idx.p1]
            founder2.df[idx.p2,] <- alt.alleles[idx.p2]
          } else if(i == length(doses)){
            if(length(idx.p1.tot)!= 0 | length(idx.p2.tot)!= 0){
              cat(length(idx.p1.tot), length(idx.p2.tot))
              for(j in 1:length(idx.p1.tot)){
                dose.idx <- sample(1:ploidy,ploidys[i])
                founder1.df[idx[idx.p1.tot][j],dose.idx] <- alt.alleles[idx[idx.p1.tot][j]]
                founder1.df[idx[idx.p1.tot][j],which(!1:ploidy == dose.idx)] <- ref.alleles[idx[idx.p1.tot][j]]
                dose.idx <- sample(1:ploidy,ploidys[i])
                founder2.df[idx[idx.p2.tot][j],dose.idx] <- alt.alleles[idx[idx.p2.tot][j]]
                founder2.df[idx[idx.p2.tot][j],which(!1:ploidy == dose.idx)] <- ref.alleles[idx[idx.p2.tot][j]]
             }
           }
          } else {
            idx.p1 <- sample(idx.p1.tot, as.numeric(size))
            idx.p2 <- sample(idx.p2.tot, as.numeric(size))
            for(j in 1:length(idx.p1)){
              dose.idx <- sample(1:ploidy,ploidys[i])
              founder1.df[idx[idx.p1][j],dose.idx] <- alt.alleles[idx[idx.p1][j]]
              founder1.df[idx[idx.p1][j],which(!1:ploidy == dose.idx)] <- ref.alleles[idx[idx.p1][j]]
              dose.idx <- sample(1:ploidy,ploidys[i])
              founder2.df[idx[idx.p2][j],dose.idx] <- alt.alleles[idx[idx.p2][j]]
              founder2.df[idx[idx.p2][j],which(!1:ploidy == dose.idx)] <- ref.alleles[idx[idx.p2][j]]
            }
            idx.p1.tot <- idx.p1.tot[-which(idx.p1.tot%in%idx.p1)]
            idx.p2.tot <- idx.p2.tot[-which(idx.p2.tot%in%idx.p2)]
          }
        }

        founder_file <- cbind(marker, founder1.df, founder2.df)
        colnames(founder_file) <- c("marker", paste0("P1_",1:ploidy), paste0("P2_",1:ploidy))

      } else if("~{cross}" == "F2"){
        founder_file <- data.frame(marker=marker,
                                   P1_1=ref_alt_alleles[,3] ,
                                   P1_2=ref_alt_alleles[,3],
                                   P2_1=ref_alt_alleles[,4],
                                   P2_2=ref_alt_alleles[,4]) # Only for diploids
      }

      write.table(founder_file, file = paste0("founders.txt"), quote=FALSE, col.names = TRUE, row.names = FALSE, sep = "\t" )

      simulated_phases <- compare_phases(founder_file)

      #create_parfile(~{seed}, 50*~{popsize})
      create_parfile(~{seed}, ~{popsize})

      create_chromfile(map_file)

     RSCRIPT

  >>>

  runtime {
    docker:"cristaniguti/reads2map:0.0.4"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "CreatePedigreeSimulatorInputs"
    mem:"~{memory_size}M"
    time:"05:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Generates the parents haplotypes when reference VCF is not provided."
  }

  output {
    File mapfile_nomap = "mapfile.txt"
    File founderfile_nomap = "founders.txt"
    File parfile_nomap = "parameters.txt"
    File chromfile_nomap = "chromosome.txt"
    File ref_alt_alleles_nomap = "ref_alt_alleles.txt"
    File simulated_phases_nomap = "simulated_phases.txt"
  }
}


# Parse pedsim output (.dat) into VCF
task ConvertPedigreeSimulationToVcf {

  input {
    Int seed
    Int depth
    File genotypes_dat
    File map_file
    File chrom_file
    File ref_alt_alleles
    Int popsize
    Int mapsize
  }

  Int disk_size = ceil(size(genotypes_dat, "GiB") * 2 + size(map_file, "GiB") + size(ref_alt_alleles, "GiB") + size(chrom_file, "GiB"))
  Int memory_size = 8000

  command <<<
    R --vanilla --no-save <<RSCRIPT

    library(onemap)
    library(Reads2MapTools)
    library(vcfR)

    mks <- read.table("~{ref_alt_alleles}", stringsAsFactors = FALSE)

    set.seed(~{seed})
    pedsim2vcf(inputfile = "~{genotypes_dat}",
               map.file = "~{map_file}",
               chrom.file = "~{chrom_file}",
               out.file = "temp.vcf",
               miss.perc = 0,
               counts = FALSE,
               pos = mks[,2],
               haplo.ref = "P1_1",
               chr = mks[,1],
               phase = TRUE,
               reference.alleles = mks[,3],
               use.as.alleles=TRUE,
              #  n_selected_loci = 1,
              #  selection_str_mean = 0.5,
              #  selection_str_var = 0.0001,
              #  pop.size = ~{popsize},
              #  selected_mks = 30,
               map.size = ~{mapsize})

    vcfR.object <- read.vcfR("temp.vcf")

    vcf_simu <- data.frame(vcfR.object@fix, vcfR.object@gt, stringsAsFactors = FALSE)

    vcf_simu[,6] <- "."
    vcf_simu[,8] <- "."

    add_head(vcf_simu, "~{seed}_~{depth}_simu.vcf")

    INDS_temp <- dimnames(vcfR.object@gt)[[2]][-1]
    inds_sele <- INDS_temp[-c(which(INDS_temp=="P1"), which(INDS_temp=="P2"))]

    progeny_dat <- vcf2progeny_haplotypes(vcfR.object = vcfR.object, ind.id = inds_sele,
                                          parent1 = "P1", parent2 = "P2",
                                          crosstype = "outcross")

    haplo_simu <- cbind(seed="~{seed}", depth="~{depth}",progeny_dat)
    vroom::vroom_write(haplo_simu, "~{seed}_~{depth}_haplo_simu.tsv.gz")

    # For RADinitio
    vcf_radinitio <- data.frame("CHROM"=1, "POS"=as.numeric(as.character(vcfR.object@fix[,2])), "ID"= ".","REF"=0, "ALT"=1,
                             "QUAL"=".", "FILTER"=".","INFO"=".",vcfR.object@gt, stringsAsFactors = FALSE)

    add_head(vcf_radinitio, "radinitio.vcf", type="radinitio")

    RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/reads2map:0.0.4"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "ConvertPedigreeSimulationToVcf"
    mem:"~{memory_size}M"
    time:"05:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Converts PedigreeSim output to VCF file. It have the option to simulate segregation distortion while this convertion is made."
  }

  output {
    File simu_vcf = "~{seed}_~{depth}_simu.vcf"
    File simu_haplo = "~{seed}_~{depth}_haplo_simu.tsv.gz"
    File radinitio_vcf = "radinitio.vcf"
  }
}


task Vcf2PedigreeSimulator{
  input {
    File vcf_file
    File? ref_map
    Int seed
    Int popsize
    String vcf_parent1
    String vcf_parent2
  }

  Int disk_size = ceil(size(vcf_file, "GiB") * 2 + size(ref_map, "GiB") + 5)
  Int memory_size = 5000

  command <<<
    R --vanilla --no-save <<RSCRIPT

    # Warning: The markers in vcf out of the reference map interval will be excluded

    library(Reads2MapTools)
    library(vcfR)

    vcf <- read.vcfR("~{vcf_file}")
    ref_map <- read.csv("~{ref_map}")

    ref_map <- remove_outlier(ref_map, thr=0) # Remove inverted markers

    # PedigreeSim inputs
    founderfile <- create_haplo(vcfR.obj = vcf, ref.map = ref_map, seed = ~{seed},
                                P1 = "~{vcf_parent1}", P2= "~{vcf_parent2}")

    ## This function generates the mapfile and the ref_alt_alleles file
    mapfile <- create_mapfile(vcf, ref_map)
    #create_parfile(~{seed}, 50*~{popsize})
    create_parfile(~{seed}, ~{popsize})
    create_chromfile(mapfile[[1]])

    ref_alt_alleles <- mapfile[[2]]
    write.table(ref_alt_alleles, file="ref_alt_alleles.txt")

    # Codifying phases for comparision with gusmap
    compare_phases(founderfile, ref_alt_alleles)

    RSCRIPT

  >>>

  runtime {
      docker:"cristaniguti/reads2map:0.0.4"
      cpu:1
      # Cloud
      memory:"~{memory_size} MiB"
      disks:"local-disk " + disk_size + " HDD"
      # Slurm
      job_name: "Vcf2PedigreeSimulator"
      mem:"~{memory_size}M"
      time:"05:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Uses information of a reference VCF to generate [PedigreeSim](https://www.wur.nl/en/show/Software-PedigreeSim.htm) inputs."
  }

  output{
      File mapfile_map = "mapfile.txt"
      File founderfile_map = "founders.txt"
      File parfile_map = "parameters.txt"
      File chromfile_map = "chromosome.txt"
      File ref_alt_alleles_map = "ref_alt_alleles.txt"
      File simulated_phases_map = "simulated_phases.txt"
  }
}

task ProduceFamiliesSeeds {
  input {
    Int number_of_families
    Int global_seed
  }

  Int disk_size = 1
  Int memory_size = 1000

  command <<<
    python <<CODE
    import random
    random.seed(~{global_seed})
    for x in range(~{number_of_families}):
        print(random.randint(1,101+x))
    CODE
  >>>

  runtime {
    docker:"python:3.7"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "ProduceFamiliesSeeds"
    mem:"~{memory_size}M"
    time:"01:00:00"
  }

  output {
    Array[Int] seeds = read_lines(stdout())
  }
}
