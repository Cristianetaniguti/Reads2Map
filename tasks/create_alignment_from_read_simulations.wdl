version 1.0

import "../structs/reads_simuS.wdl"
import "../structs/alignment_struct.wdl"
import "alignment.wdl" as alg

workflow CreateAlignmentFromSimulation {
    input {
        ReferenceFasta references
        Family family
        Profiles profiles
    }


  call GenerateAlternativeGenome {
    input:
      seed       = family.seed,
      ref_genome = references.ref_fasta
  }

  call CreatePedigreeSimulatorInputs {
    input:
      seed       = family.seed,
      snps       = GenerateAlternativeGenome.snps,
      indels     = GenerateAlternativeGenome.indels,
      cmBymb     = family.cmBymb,
      ref        = references.ref_fasta,
      ref_fai    = references.ref_fasta_index,
      cross      = family.cross,
      popsize    = family.popsize,
      ploidy     = family.ploidy,
      doses      = family.doses
  }

  call RunPedigreeSimulator {
    input:
      mapfile     = CreatePedigreeSimulatorInputs.mapfile,
      founderfile = CreatePedigreeSimulatorInputs.founderfile,
      chromfile   = CreatePedigreeSimulatorInputs.chromfile,
      parfile     = CreatePedigreeSimulatorInputs.parfile
  }

  call ConvertPedigreeSimulationToVcf {
    input:
      seed       = family.seed,
      depth      = family.depth,
      genotypes_dat = RunPedigreeSimulator.genotypes_dat,
      map_file      = CreatePedigreeSimulatorInputs.mapfile,
      chrom_file    = CreatePedigreeSimulatorInputs.chromfile,
      tot_mks       = CreatePedigreeSimulatorInputs.tot_mks
  }

  call GenerateSampleNames {
    input:
      simulated_vcf = ConvertPedigreeSimulationToVcf.simu_vcf
  }

  scatter (sampleName in GenerateSampleNames.names) {

    call RunVcf2diploid {
      input:
        sampleName = sampleName,
        ref_genome = references.ref_fasta,
        simu_vcf   = ConvertPedigreeSimulationToVcf.simu_vcf
    }

    call SimulateRADseq {
      input:
        enzyme           = family.enzyme,
        sampleName       = sampleName,
        maternal_genomes = RunVcf2diploid.maternal_genomes,
        paternal_genomes = RunVcf2diploid.paternal_genomes
    }

    call SimulateIlluminaReads {
      input:
        maternal_trim = SimulateRADseq.maternal_trim,
        paternal_trim = SimulateRADseq.paternal_trim,
        sampleName    = sampleName,
        depth         = family.depth,
        base_calling  = profiles.base_calling,
        indel_error   = profiles.indel_error,
        gc_bias       = profiles.gc_bias
    }

    call alg.RunBwaAlignment {
      input:
        sampleName = sampleName,
        reads1     = [SimulateIlluminaReads.reads1],
        libraries  = ["artificial"],
        ref        = references.ref_fasta,
        geno_amb   = references.ref_amb,
        geno_ann   = references.ref_ann,
        geno_bwt   = references.ref_bwt,
        geno_pac   = references.ref_pac,
        geno_sa    = references.ref_sa
    }
  }

  output {
      Array[Alignment] alignments = RunBwaAlignment.algn
      Array[File] bam = RunBwaAlignment.bam
      Array[File] bai = RunBwaAlignment.bai
      File total_markers = CreatePedigreeSimulatorInputs.tot_mks
      Array[File] maternal_trim = SimulateRADseq.maternal_trim
      Array[String] names = GenerateSampleNames.names
      File true_vcf = ConvertPedigreeSimulationToVcf.simu_vcf
      File simu_haplo = ConvertPedigreeSimulationToVcf.simu_haplo
      File real_phases = CreatePedigreeSimulatorInputs.real_phases
  }

}


# Creates homologous genome with some variation
# specified with -s and -d
task GenerateAlternativeGenome {
  input {
    Int seed
    File ref_genome
  }

  command <<<
    /pirs/src/pirs/pirs diploid ~{ref_genome} -s 0.001 -d 0.0001 -v 0 -o alt --random-seed ~{seed}

  >>>

  runtime {
    docker: "taniguti/pirs-ddrad-cutadapt"
    # mem:"--nodes=1"
    cpu:1
    time:"48:00:00"
  }

  output {
    File alt_fasta = "alt.snp.indel.fa"
    File indels = "alt.indel.lst"
    File snps = "alt.snp.lst"
  }
}

task CreatePedigreeSimulatorInputs {
  input {
    File snps
    File indels
    Float cmBymb
    File ref
    File ref_fai
    Int seed
    Int popsize
    Int ploidy # Only ploidy 2 available
    File doses
    String cross
  }


  command <<<

      R --vanilla --no-save <<RSCRIPT

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
      command  <- c(paste("samtools faidx ~{ref}"),paste0(indels[1,1],":",pos.ref,"-",pos.ref))
      bases <- system(paste0(command, collapse = " "), intern = T)

      bases.bf <- matrix(bases, ncol=2, byrow = T)[,2]
      alt <- bases.bf
      tmp <- paste0(bases.bf[which(sinal=="+")], indels[,6][which(sinal=="+")])
      alt[which(sinal=="+")] <- tmp
      ref <- bases.bf
      tmp <- paste0(bases.bf[which(sinal=="-")], indels[,6][which(sinal=="-")])
      ref[which(sinal=="-")] <- tmp

      # the position in the vcf and in the map are according with the reference genome
      tot.mks <- data.frame(chr = c(snps[,1], indels[,1]), pos = c(snps[,2], pos.ref),
                            ref = c(snps[,4], ref), alt = c(snps[,5],alt), stringsAsFactors = F)


      tot.mks <- tot.mks[order(tot.mks[,2]),]
      write.table(tot.mks, file="tot_mks.txt")
      n.marker <- dim(tot.mks)[1]

      ## Map file
      # Marker names
      marker1 <- "M"
      marker2 <- 1:n.marker
      marker2 <- sprintf("%03d", marker2)
      marker <-paste0(marker1,marker2)

      # Chromossome and position
      pos.map <- (tot.mks[,2]/1000000) * ~{cmBymb}
      map_file <- data.frame(marker=marker, chromosome=tot.mks[,1], position= pos.map)
      write.table(map_file, file = paste0("mapfile.txt"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

      ## Founderfile
      ref.alleles <- tot.mks[,3]
      alt.alleles  <- tot.mks[,4]

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
        founder_file <- data.frame(marker=marker, P1_1=tot.mks[,3] , P1_2=tot.mks[,3], P2_1=tot.mks[,4], P2_2=tot.mks[,4]) # Only for diploids
      }

      write.table(founder_file, file = paste0("founders.txt"), quote=FALSE, col.names = TRUE, row.names = FALSE, sep = "\t" )

      ## Real Phases for comparisions with Gusmap
      founder <- founder_file[,-1]
      real_phases <- rep(NA, dim(founder)[1])
      real_phases[which(founder[,1] == founder[,3] & founder[,2] == founder[,4])] <- 17 # 1 and 4
      real_phases[which(founder[,1] == founder[,4] & founder[,2] == founder[,3])] <- 18 # 2 and 3
      real_phases[which(founder[,1] == founder[,3] & founder[,1] == founder[,4] & founder[,1] != founder[,2])] <- 19 # 5 and 8
      real_phases[which(founder[,2] == founder[,3] & founder[,2] == founder[,4] & founder[,1] != founder[,2])] <- 20 # 6 and 7
      real_phases[which(founder[,1] == founder[,2] & founder[,1] == founder[,3] & founder[,1] != founder[,4])] <- 21 # 9 and 12
      real_phases[which(founder[,1] == founder[,2] & founder[,1] == founder[,4] & founder[,1] != founder[,3])] <- 22 # 10 and 11
      real_phases[which(founder[,1] == founder[,2] & founder[,1] == founder[,3] & founder[,1] == founder[,4])] <- 23 # 13 and 16
      real_phases[which(founder[,1] == founder[,2] & founder[,3] == founder[,4] & founder[,1] != founder[,3])] <- 24 # 14 and 15

      real_phases <- data.frame(pos=tot.mks[,2], real_phases)

      write.table(real_phases, file = paste0("real_phases.txt"))

      ## Parameters file
      parameter <- paste0("PLOIDY = ~{ploidy}
                           MAPFUNCTION = HALDANE
                           MISSING = NA
                           CHROMFILE = chromosome.txt
                           POPTYPE = ~{cross}
                           SEED = ~{seed}
                           POPSIZE = ~{popsize}
                           MAPFILE = mapfile.txt
                           FOUNDERFILE = founderfile.txt
                           OUTPUT = sim")

      write.table(parameter, file = paste0("parameters.txt"), quote=FALSE, col.names = FALSE, row.names = FALSE, sep = "\t" )
      chrom <- data.frame("chromosome"= "Chr10", "length"= pos.map[which.max(pos.map)], "centromere"=pos.map[which.max(pos.map)]/2, "prefPairing"= 0.0, "quadrivalents"=0.0)
      write.table(chrom, file= "chromosome.txt", quote = F, col.names = T, row.names = F, sep= "\t")
     RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/r-samtools"
    # mem:"--nodes=1"
    cpu:1
    time:"24:00:00"
  }

  output {
    File mapfile = "mapfile.txt"
    File founderfile = "founders.txt"
    File parfile = "parameters.txt"
    File chromfile = "chromosome.txt"
    File tot_mks = "tot_mks.txt"
    File real_phases = "real_phases.txt"
  }

}

task RunPedigreeSimulator {
  input {
    File mapfile
    File founderfile
    File chromfile
    File parfile
  }

  command <<<
    set -e
    sed -i 's+chromosome.txt+~{chromfile}+g' ~{parfile}
    sed -i 's+mapfile.txt+~{mapfile}+g' ~{parfile}
    sed -i 's+founderfile.txt+~{founderfile}+g' ~{parfile}
    java -jar /usr/jars/PedigreeSim.jar ~{parfile}

  >>>

  runtime {
    docker: "taniguti/java-in-the-cloud"
    # mem:"--nodes=1"
    cpu:1
    time:"24:00:00"
  }

  output {
    File genotypes_dat = "sim_genotypes.dat"
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
    File tot_mks
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT

    library(onemap)
    library(vcfR)

    mks <- read.table("~{tot_mks}", stringsAsFactors = FALSE)
    pos <- mks[,2]
    chr <- mks[,1]

    pedsim2vcf(inputfile = "~{genotypes_dat}",
      map.file = "~{map_file}",
      chrom.file = "~{chrom_file}",
      out.file = "~{seed}_~{depth}_simu.vcf",
      miss.perc = 0, counts = FALSE,pos = pos, haplo.ref = "P1_1",
      chr = chr, phase = TRUE)

    vcfR.object <- read.vcfR("~{seed}_~{depth}_simu.vcf")
    INDS_temp <- dimnames(vcfR.object@gt)[[2]][-1]
    inds_sele <- INDS_temp[-c(which(INDS_temp=="P1"), which(INDS_temp=="P2"))]

    progeny_dat <- vcf2progeny_haplotypes(vcfR.object = vcfR.object, ind.id = inds_sele,
                                          parent1 = "P1", parent2 = "P2",
                                          crosstype = "outcross")

    haplo_simu <- cbind(seed="~{seed}", depth="~{depth}",progeny_dat)
    saveRDS(haplo_simu, file = "~{seed}_~{depth}_haplo_simu.rds")

    RSCRIPT

  >>>

  runtime {
    docker: "gcr.io/taniguti-backups/onemap:v1"
    # mem:"--nodes=1"
    cpu:1
    time:"48:00:00"
  }

  output {
    File simu_vcf = "~{seed}_~{depth}_simu.vcf"
    File simu_haplo = "~{seed}_~{depth}_haplo_simu.rds"
  }
}

# Insert into a fasta sequence the variants present in a VCF file
task RunVcf2diploid {
  input {
    String sampleName
    File ref_genome
    File simu_vcf
  }

  command <<<
    java -jar /usr/jars/vcf2diploid.jar -id ~{sampleName} -chr ~{ref_genome} -vcf ~{simu_vcf}

  >>>

  runtime {
    docker: "taniguti/java-in-the-cloud"
    # mem:"--nodes=1"
    cpu:1
    time:"24:00:00"
  }

  output {
    File maternal_genomes = "Chr10_${sampleName}_maternal.fa"
    File paternal_genomes = "Chr10_${sampleName}_paternal.fa"
  }

}

# It will always produce P1, P2, F1 and then F2_00X, where
# X will increase from 1 to samples
task GenerateSampleNames {

  input {
    File simulated_vcf
  }

  command <<<
    export PATH=$PATH:/opt/conda/bin

    python <<CODE
    from pysam import VariantFile

    bcf_in = VariantFile("~{simulated_vcf}")

    for i in bcf_in.header.samples:
        print(i)
    CODE

  >>>

  runtime {
    docker: "taniguti/miniconda-alpine"
    # mem:"--mem-per-cpu=24042"
    cpu:1
    time:"24:00:00"
  }

  output {
    Array[String] names = read_lines(stdout())
  }
}


# Simulates RADseq experiment where certain enzyme is
# used to fragment the sequence and then the first X
# bases of each resulting fragment is sequenced.
task SimulateRADseq {
  input {
    String enzyme
    String sampleName
    File maternal_genomes
    File paternal_genomes
  }

  command <<<
    set -e
    /ddRADseqTools/Package/rsitesearch.py \
      --genfile=~{maternal_genomes} \
      --fragsfile=~{sampleName}_maternal_fragments.fasta \
      --rsfile=/ddRADseqTools/Package/restrictionsites.txt \
      --enzyme1=~{enzyme} \
      --enzyme2=~{enzyme} \
      --minfragsize=202 \
      --maxfragsize=500 \
      --fragstfile=~{sampleName}_maternal_statistics.txt \
      --fragstinterval=25 \
      --plot=NO \
      --verbose=YES \
      --trace=NO

    /ddRADseqTools/Package/rsitesearch.py \
      --genfile=~{paternal_genomes} \
      --fragsfile=~{sampleName}_paternal_fragments.fasta \
      --rsfile=/ddRADseqTools/Package/restrictionsites.txt \
      --enzyme1=~{enzyme} \
      --enzyme2=~{enzyme} \
      --minfragsize=202 \
      --maxfragsize=500 \
      --fragstfile=~{sampleName}_paternal_statistics.txt \
      --fragstinterval=25 \
      --plot=NO \
      --verbose=YES \
      --trace=NO

    cutadapt -l 202 \
      -o ~{sampleName}_maternal_trim.fa \
      ~{sampleName}_maternal_fragments.fasta

    cutadapt -l 202 \
      -o ~{sampleName}_paternal_trim.fa \
      ~{sampleName}_paternal_fragments.fasta

  >>>

  runtime {
    docker: "taniguti/pirs-ddrad-cutadapt"
    # mem:"--nodes=1"
    cpu:1
    time:"48:00:00"
  }

  output {
    File maternal_frags = "${sampleName}_maternal_fragments.fasta"
    File paternal_frags = "${sampleName}_paternal_fragments.fasta"
    File maternal_stats = "${sampleName}_maternal_statistics.txt"
    File paternal_stats = "${sampleName}_paternal_statistics.txt"
    File maternal_trim = "${sampleName}_maternal_trim.fa"
    File paternal_trim = "${sampleName}_paternal_trim.fa"
  }
}

task SimulateIlluminaReads {

  input {
    File maternal_trim
    File paternal_trim
    Int depth
    String sampleName
    File base_calling
    File indel_error
    File gc_bias
  }

  command <<<
    set -e
    /pirs/src/pirs/pirs simulate \
      --diploid ~{maternal_trim} ~{paternal_trim} \
      --read-len=100 \
      --coverage=~{depth} \
      --insert-len-mean=150 \
      --output-prefix=~{sampleName} \
      --output-file-type=gzip \
      --threads=20 \
      --base-calling-profile=~{base_calling} \
      --indel-error-profile=~{indel_error} \
      --gc-bias-profile=~{gc_bias}

  >>>

  runtime {
    docker: "taniguti/pirs-ddrad-cutadapt"
    maxRetries: 5
    # mem:"--nodes=1"
    time:"48:00:00"
    cpu:20
  }

  output {
    File reads1 = "${sampleName}_100_150_1.fq.gz"
    File reads2 = "${sampleName}_100_150_2.fq.gz"
  }
}
