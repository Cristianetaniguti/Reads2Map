version 1.0

import "../structs/reads_simuS.wdl"
import "../structs/alignment_struct.wdl"
import "alignment.wdl" as alg

workflow CreateAlignmentFromSimulation {
    input {
        Reference references
        Family family
        Sequencing sequencing
    }

  if (!defined(sequencing.vcf)){
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
  }

  if (defined(sequencing.vcf)){
    call Vcf2PedigreeSimulator {
      input:
        vcf_file = sequencing.vcf,
        ref_map = sequencing.ref_map,
        parent1 = family.parent1,
        parent2 = family.parent2,
        seed = family.seed,
        popsize = family.popsize,
        cmBymb = family.cmBymb
    }
  }

  File mapfile = select_first([Vcf2PedigreeSimulator.mapfile, CreatePedigreeSimulatorInputs.mapfile])
  File founderfile = select_first([Vcf2PedigreeSimulator.founderfile, CreatePedigreeSimulatorInputs.founderfile])
  File chromfile = select_first([Vcf2PedigreeSimulator.chromfile, CreatePedigreeSimulatorInputs.chromfile])
  File parfile = select_first([Vcf2PedigreeSimulator.parfile, CreatePedigreeSimulatorInputs.parfile])

  call RunPedigreeSimulator {
    input:
      mapfile     = mapfile,
      founderfile = founderfile,
      chromfile   = chromfile,
      parfile     = parfile
  }

  call ConvertPedigreeSimulationToVcf {
    input:
      seed       = family.seed,
      depth      = family.depth,
      genotypes_dat = RunPedigreeSimulator.genotypes_dat,
      map_file      = CreatePedigreeSimulatorInputs.mapfile,
      chrom_file    = CreatePedigreeSimulatorInputs.chromfile,
      ref_alt_alleles       = CreatePedigreeSimulatorInputs.ref_alt_alleles
  }

  call GenerateSampleNames {
    input:
      simulated_vcf = ConvertPedigreeSimulationToVcf.simu_vcf
  }

  if(sequencing.type == "WGS" || sequencing.type == "exome"){

    call SimuscopProfile {
      input:
        type = sequecing.type,
        emp_bam = sequencing.emp_bam,
        vcf     = ConvertPedigreeSimulationToVcf.simu_vcf,
        references = references
    }

    scatter (sampleName in GenerateSampleNames.names) {

      call SimuscopSimulation {
        input:
          type          = sequencing.type,
          sampleName    = sampleName,
          depth         = family.depth,
          emp_bam       = sequencing.emp_bam,
          vcf           = ConvertPedigreeSimulationToVcf.simu_vcf,
          references    = references,
          chrom         = sequencing.chromosome,
          profile       = SimuscopProfile.profile
      }
    }
  }

  # Two option of RADseq
  if(sequencing.type == "sdRAD" || sequencing.type == "ddRAD"){
    call RADinitioSimulation{
      input:
        depth          = sequencing.depth,
        enzyme1        = sequencing.enzyme1,
        enzyme2        = sequencing.enzyme2,
        vcf_simu       = ConvertPedigreeSimulationToVcf.simu_vcf,
        references     = references,
        pcr_cycles     = sequencing.pcr_cycles,
        insert_size    = sequencing.insert_size,   
        isert_size_dev = sequencing.isert_size_dev,
        read_length    = sequencing.read_length,
        library_type   = sequencing.type,
        sampleName     = sampleName
    }
  }

  Array[File] fastqs = select_first([RADinitioSimulation.fastq, SimuscopSimulation.fastq])

  scatter (sampleName in GenerateSampleNames.names) {
    call alg.RunBwaAlignment {
      input:
        sampleName = sampleName,
        reads1     = fastqs,
        libraries  = ["artificial"],
        references = references
    }
  }

  output {
      Array[Alignment] alignments = RunBwaAlignment.algn
      Array[File] bam = RunBwaAlignment.bam
      Array[File] bai = RunBwaAlignment.bai
      File total_markers = CreatePedigreeSimulatorInputs.ref_alt_alleles
      Array[File] maternal_trim = SimulateRADseq.maternal_trim
      Array[String] names = GenerateSampleNames.names
      File true_vcf = ConvertPedigreeSimulationToVcf.simu_vcf
      File simu_haplo = ConvertPedigreeSimulationToVcf.simu_haplo
      File simulated_phases = CreatePedigreeSimulatorInputs.simulated_phases
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
    /pirs/src/pirs/pirs diploid ~{ref_genome} -s 0.0133 -d 0.0022 -v 0 -o alt --random-seed ~{seed}

  >>>

  runtime {
    docker: "taniguti/pirs-ddrad-cutadapt"
    mem:"--nodes=1"
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
      write.table(ref_alt_alleles, file="ref_alt_alleles.txt")
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
        founder_file <- data.frame(marker=marker, P1_1=ref_alt_alleles[,3] , P1_2=ref_alt_alleles[,3], P2_1=ref_alt_alleles[,4], P2_2=ref_alt_alleles[,4]) # Only for diploids
      }

      write.table(founder_file, file = paste0("founders.txt"), quote=FALSE, col.names = TRUE, row.names = FALSE, sep = "\t" )

      source("/opt/scripts/vcf2pedigreeSim.R")
      simulated_phases <- compare_phases(founder_file)

      create_parfile(~{seed}, ~{popsize})

      create_chromfile(map_file)

     RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/r-samtools"
    mem:"--nodes=1"
    cpu:1
    time:"24:00:00"
  }

  output {
    File mapfile = "mapfile.txt"
    File founderfile = "founders.txt"
    File parfile = "parameters.txt"
    File chromfile = "chromosome.txt"
    File ref_alt_alleles = "ref_alt_alleles.txt"
    File simulated_phases = "simulated_phases.txt"
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
    mem:"--nodes=1"
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
    File ref_alt_alleles
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT

    library(onemap)
    library(vcfR)

    mks <- read.table("~{ref_alt_alleles}", stringsAsFactors = FALSE)
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
    docker: "cristaniguti/onemap_workflows"
    mem:"--nodes=1"
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
    mem:"--nodes=1"
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
    mem:"--mem-per-cpu=24042"
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
# Deprecated!!!
task SimulateRADseq {
  input {
    String enzyme1
    String enzyme2
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
      --enzyme1=~{enzyme1} \
      --enzyme2=~{enzyme2} \
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
      --enzyme1=~{enzyme1} \
      --enzyme2=~{enzyme2} \
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
    mem:"--nodes=1"
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

# Deprecated
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
    mem:"--nodes=1"
    time:"48:00:00"
    cpu:20
  }

  output {
    File reads1 = "${sampleName}_100_150_1.fq.gz"
    File reads2 = "${sampleName}_100_150_2.fq.gz"
  }
}


task Vcf2PedigreeSimulator{
  input{
    File vcf_file
    File? ref_map
    String parent1
    String parent2
    Int seed
    Int popsize
    Float? cmBymb
  }

    command <<<
    R --vanilla --no-save <<RSCRIPT

    source("/opt/scripts/vcf2pedigreeSim.R")

    vcf <- read.vcfR("~{vcf_file}")
    ref_map <- read.csv("~{ref_map}")

    # PedigreeSim inputs
    founderfile <- create_haplo(vcfR.obj = vcf, seed = 1010, 
                                P1 = "~{parent1}", P2= "~{parent2}")

    mapfile <- create_mapfile(vcf, ref_map)
    create_parfile(~{seed}, ~{popsize})
    create_chromfile(mapfile)

    # Store ref and alt alleles in file
    ref_alt_alleles <- data.frame(chr = vcf@fix[,1], 
                                  pos = vcf@fix[,2], 
                                  ref = vcf@fix[,4], 
                                  alt = vcf@fix[,5], stringsAsFactors = F)

    write.table(ref_alt_alleles, file="ref_alt_alleles.txt")

    # Codifying phases for comparision with gusmap
    compare_phases(founderfile)


    RSCRIPT

    runtime {
      docker: "cristaniguti/r-samtools"
      mem:"--nodes=1"
      cpu:1
      time:"24:00:00"
     }

    output{
      File mapfile = "mapfile.txt"
      File founderfile = "founders.txt"
      File parfile = "parameters.txt"
      File chromfile = "chromosome.txt"
      File ref_alt_alleles = "ref_alt_alleles.txt"
      File simulated_phases = "simulated_phases.txt"
    }

}

task SimuscopProfile{
  inputs{
    String type          
    File   emp_bam     
    File   vcf       
    File   references
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT

    library(vcfR)
    library(simuscopR)

    if("~{type}" == "exome"){

      system(paste0("bamtobed -i ~{emp_bam} > bed_file"))

      seqToProfile("~{emp_bam}", "bed_file", "~{vcf}",
             "~{references.ref_fasta}", "profile")
      
    } else {
      seqToProfile("~{emp_bam}", vcf.file =  "~{vcf}",
             reference = "~{references.ref_fasta}",  out.profile = "sample.profile")
    }

    RSCRIPT

    runtime {
      docker: "cristaniguti/simuscopr"
      mem:"--nodes=1"
      cpu:1
      time:"24:00:00"
    }

    output{
      File profile = "sample.profile"
    }

}

task SimuscopSimulation{
 inputs{
    String type          
    String sampleName
    Int    depth    
    File   emp_bam     
    File   vcf       
    File   references
    String chrom
    File   profile
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
    vcfR.object <- read.vcfR("~{vcf}")

    variants <- vcf2variants(vcfR.object, sample = "~{sampleName}", chrom = "~{chrom}")

    write.table(variants$SNVs, file = "SNVs.txt", sep = "\t", quote = F, col.names = F, row.names = F)
    write.table(variants$indels, file = "indels.txt", sep = "\t", quote = F, col.names = F, row.names = F)
    write.table(variants$insertions, file = "insertions.txt", sep = "\t", quote = F, col.names = F, row.names = F)

    system("cat SNVs.txt indels.txt insertions.txt > variants.txt")

    if("~{type}" == "exome"){
      simuReads(ref = "~{references.ref_fasta}",
            profile = "~{profile}",
            variation = "variants.txt",
            target = "bed_file",
            name = "~{sanmpleName}",
            output = ".",
            layout = "SE",
            threads = 6,
            verbose = 1,
            coverage = ~{depth})
    } else {
      simuReads(ref = "~{references.ref_fasta}",
            profile = "profile",
            variation = "variants.txt",
            name = "~{sanmpleName_ref}",
            output = ".",
            layout = "SE", # only single-end by now
            threads = 6,
            verbose = 1,
            coverage = ~{depth})
    }
    

  RSCRIPT

  runtime {
    docker: "cristaniguti/simuscopr"
    mem:"--nodes=1"
    cpu:1
    time:"24:00:00"
   }

  output {
    File fastq = "~{sampleName}.fq"
  }
}


task RADinitioSimulation{
  inputs{
    File vcf_simu
    String enzyme1
    String enzyme2
    File references
    Int depth
    Int? insert_size
    Int? isert_size_dev
    Int? pcr_cycles
    Int? read_length
    String library_type
    String sampleName
  }

  command <<<
    cat ~{references.ref_fasta_index} | cut -f 1 | head -n 10 | tail -n 1 > ./chrom.list

    # Makes for individual sample for parallelization
    echo -e ~{sampleName}'\t'pop0 > popmap.tsv

    vcftools --vcf ~{simu_vcf} --indv ~{sampleName} --recode --out ~{sampleName}

    mkdir msprime_vcfs ref_loci_vars
    gzip ~{sampleName}.recode.vcf 
    cp ~{sampleName}.recode.vcf.gz msprime_vcfs/~{chrom}.vcf.gz
    mv ~{sampleName}.recode.vcf.gz ref_loci_vars/ri_master.vcf.gz

    mkdir simu_inputs results
    mv msprime_vcfs ref_loci_vars popmap.tsv simu_inputs

    radinitio --make-library-seq \
              --genome ~{references.ref_fasta} \
              --chromosomes chrom.list \
              --out-dir results/ \
              --make-pop-sim-dir simu_inputs/ \
              --library-type ~{library_type} \
              --enz ~{enzyme1} \
              ~{"--enz2 " + enzyme2} \
              --insert-mean ~{default="350" insert_size} \
              --insert-stdev ~{default="35" insert_size_dev} \
              --pcr-cycles ~{default="9" pcr_cycles} \
              --coverage ~{default="20" depth} \
              --read-length ~{default="150" read_length}
    
    # Add fake phred score of 40 (H in Illumina 1.8+ Phred+33)
    seqtk seq -F 'H' results/rad_reads/~{sampleName}.1.fq.gz > ~{sampleName}.fq

  >>>

  runtime{
    docker: "cristaniguti/radioinitio"
    mem:"--nodes=1"
    cpu:1
    time:"24:00:00"
  }

  output{
    File fastq = "~{sampleName}.fq"
  }
}
