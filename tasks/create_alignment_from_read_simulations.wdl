version 1.0

import "../structs/reads_simuS.wdl"
import "../structs/alignment_struct.wdl"
import "./alignment.wdl" as alg

workflow CreateAlignmentFromSimulation {
    input {
        Reference references
        Family family
        Sequencing sequencing
        Int max_cores
    }

  # User can provide specific variants in VCF file
  # If not, use pirs to simulate
  if (!defined(sequencing.emp_vcf)){
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

  if (defined(sequencing.emp_vcf)){
    call Vcf2PedigreeSimulator {
      input:
        vcf_file = sequencing.emp_vcf,
        ref_map = sequencing.ref_map,
        seed = family.seed,
        popsize = family.popsize,
        vcf_parent1 = sequencing.vcf_parent1,
        vcf_parent2 = sequencing.vcf_parent2
    }
  }

  File mapfile_sele = select_first([Vcf2PedigreeSimulator.mapfile_map, CreatePedigreeSimulatorInputs.mapfile_nomap])
  File founderfile_sele = select_first([Vcf2PedigreeSimulator.founderfile_map, CreatePedigreeSimulatorInputs.founderfile_nomap])
  File chromfile_sele = select_first([Vcf2PedigreeSimulator.chromfile_map, CreatePedigreeSimulatorInputs.chromfile_nomap])
  File parfile_sele = select_first([Vcf2PedigreeSimulator.parfile_map, CreatePedigreeSimulatorInputs.parfile_nomap])
  File ref_alt_alleles_sele = select_first([Vcf2PedigreeSimulator.ref_alt_alleles_map, CreatePedigreeSimulatorInputs.ref_alt_alleles_nomap])
  File simulated_phases_sele = select_first([Vcf2PedigreeSimulator.simulated_phases_map, CreatePedigreeSimulatorInputs.simulated_phases_nomap])


  call RunPedigreeSimulator {
    input:
      mapfile     = mapfile_sele,
      founderfile = founderfile_sele,
      chromfile   = chromfile_sele,
      parfile     = parfile_sele
  }

  #
  call ConvertPedigreeSimulationToVcf {
    input:
      seed            = family.seed,
      depth           = sequencing.depth,
      genotypes_dat   = RunPedigreeSimulator.genotypes_dat,
      map_file        = mapfile_sele,
      chrom_file      = chromfile_sele,
      ref_alt_alleles = ref_alt_alleles_sele
  }

  call GenerateSampleNames {
    input:
      simulated_vcf = ConvertPedigreeSimulationToVcf.simu_vcf
  }

  if(sequencing.library_type == "WGS" || sequencing.library_type == "exome"){

    call SimuscopProfile {
      input:
        library_type = sequencing.library_type,
        emp_bam = sequencing.emp_bam,
        vcf     = ConvertPedigreeSimulationToVcf.simu_vcf,
        references = references
    }

    scatter (sampleName in GenerateSampleNames.names) {

      call SimuscopSimulation {
        input:
          library_type  = sequencing.library_type,
          sampleName    = sampleName,
          depth         = sequencing.depth,
          emp_bam       = sequencing.emp_bam,
          vcf           = ConvertPedigreeSimulationToVcf.simu_vcf,
          references    = references,
          chrom         = sequencing.chromosome,
          profile       = SimuscopProfile.profile
      }
    }
  }


    # Two option of RADseq
    # The samples need to be simulated together, otherwise they will be all heterozygous
    if(sequencing.library_type == "sdRAD" || sequencing.library_type == "ddRAD"){
      call RADinitioSimulation{
        input:
          depth          = sequencing.depth,
          depth_parents  = sequencing.depth_parents,
          enzyme1        = sequencing.enzyme1,
          enzyme2        = sequencing.enzyme2,
          simu_vcf       = ConvertPedigreeSimulationToVcf.simu_vcf,
          radinitio_vcf  = ConvertPedigreeSimulationToVcf.radinitio_vcf,
          references     = references,
          pcr_cycles     = sequencing.pcr_cycles,
          insert_size    = sequencing.insert_size,
          insert_size_dev = sequencing.insert_size_dev,
          read_length    = sequencing.read_length,
          library_type   = sequencing.library_type,
          chrom          = sequencing.chromosome,
          names          = GenerateSampleNames.names
      }
    }


  Array[File] fastq = select_first([RADinitioSimulation.fastq_rad, SimuscopSimulation.fastq_seq])

  call alg.RunBwaAlignmentSimu {
    input:
      reads     = fastq,
      references = references,
      max_cores = max_cores,
      rm_dupli = sequencing.rm_dupli
  }

  output {
      Array[File] bam = RunBwaAlignmentSimu.bam
      Array[File] bai = RunBwaAlignmentSimu.bai
      File ref_alt_alleles = ref_alt_alleles_sele
      Array[String] names = GenerateSampleNames.names
      File true_vcf = ConvertPedigreeSimulationToVcf.simu_vcf
      File simu_haplo = ConvertPedigreeSimulationToVcf.simu_haplo
      File simulated_phases = simulated_phases_sele
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
    docker: "cristaniguti/pirs-ddrad-cutadapt:0.0.1"
    memory: "4 GB"
    cpu: 1
    preemptible: 3
    disks: "local-disk " + 5 + " HDD"
  }

  output {
    File alt_fasta = "alt.snp.indel.fa"
    File indels = "alt.indel.lst"
    File snps = "alt.snp.lst"
  }
}

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

  Int disk_size = ceil(size(ref, "GB") + size(snps, "GB") + size(indels, "GB") + 5)

  command <<<

      R --vanilla --no-save <<RSCRIPT
      library(onemapUTILs)
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

      create_parfile(~{seed}, ~{popsize})

      create_chromfile(map_file)

     RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/reads2map:0.0.1"
    memory:"8 GB"
    cpu: 2
    preemptible: 3
    disks: "local-disk " + disk_size + " HDD"
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

task RunPedigreeSimulator {
  input {
    File mapfile
    File founderfile
    File chromfile
    File parfile
  }

  Int disk_size = ceil(size(mapfile, "GB") + size(founderfile, "GB") + 5)

  command <<<
    set -e
    sed -i 's+chromosome.txt+~{chromfile}+g' ~{parfile}
    sed -i 's+mapfile.txt+~{mapfile}+g' ~{parfile}
    sed -i 's+founderfile.txt+~{founderfile}+g' ~{parfile}
    java -jar /usr/jars/PedigreeSim.jar ~{parfile}

  >>>

  runtime {
    docker: "cristaniguti/java-in-the-cloud:0.0.1"
    memory: "3 GB"
    cpu: 1
    preemptible: 3
    disks: "local-disk " + disk_size + " HDD"
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

  Int disk_size = ceil(size(genotypes_dat, "GB") + size(map_file, "GB") + 5 )

  command <<<
    R --vanilla --no-save <<RSCRIPT

    library(onemap)
    library(onemapUTILS)
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
               reference.alleles = mks[,3]
               use.as.alleles=TRUE)

    vcfR.object <- read.vcfR("temp.vcf")

    change_header(vcfR.object, "~{seed}_~{depth}_simu.vcf")

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
    docker: "cristaniguti/reads2map:0.0.1"
    memory: "4 GB"
    cpu:1
    preemptible: 3
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File simu_vcf = "~{seed}_~{depth}_simu.vcf"
    File simu_haplo = "~{seed}_~{depth}_haplo_simu.tsv.gz"
    File radinitio_vcf = "radinitio.vcf"
  }
}

# Insert into a fasta sequence the variants present in a VCF file
task RunVcf2diploid {
  input {
    String sampleName
    File ref_genome
    File simu_vcf
  }

  Int disk_size = ceil(size(ref_genome, "GB") + size(simu_vcf, "GB") + 2)

  command <<<
    java -jar /usr/jars/vcf2diploid.jar -id ~{sampleName} -chr ~{ref_genome} -vcf ~{simu_vcf}

  >>>

  runtime {
    docker: "cristaniguti/java-in-the-cloud:0.0.1"
    memory: "3 GB"
    cpu: 1
    preemptible: 3
    disks: "local-disk " + disk_size + " HDD"
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

  Int disk_size = ceil(size(simulated_vcf, "GB") + 2)

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
    docker: "cristaniguti/miniconda-alpine:0.0.1"
    memory: "1 GB"
    cpu: 1
    preemptible: 3
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    Array[String] names = read_lines(stdout())
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

  Int disk_size = ceil(size(vcf_file, "GB") +  5)

  command <<<
    R --vanilla --no-save <<RSCRIPT

    # Warning: The markers in vcf out of the reference map interval will be excluded

    library(onemapUTILS)
    library(vcfR)

    vcf <- read.vcfR("~{vcf_file}")
    ref_map <- read.csv("~{ref_map}")

    ref_map <- remove_outlier(ref_map, thr=0) # Remove inverted markers

    # PedigreeSim inputs
    founderfile <- create_haplo(vcfR.obj = vcf, ref.map = ref_map, seed = ~{seed},
                                P1 = "~{vcf_parent1}", P2= "~{vcf_parent2}")

    ## This function generates the mapfile and the ref_alt_alleles file
    mapfile <- create_mapfile(vcf, ref_map)
    create_parfile(~{seed}, ~{popsize})
    create_chromfile(mapfile[[1]])

    ref_alt_alleles <- mapfile[[2]]
    write.table(ref_alt_alleles, file="ref_alt_alleles.txt")

    # Codifying phases for comparision with gusmap
    compare_phases(founderfile, ref_alt_alleles)

    RSCRIPT

  >>>

  runtime {
      docker: "cristaniguti/reads2map:0.0.1"
      memory: "4 GB"
      cpu:1
      preemptible: 3
      disks: "local-disk " + disk_size + " HDD"
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

task SimuscopProfile{
  input {
    String library_type
    File?  emp_bam
    File   vcf
    Reference  references
  }

  Int disk_size = ceil(size(vcf, "GB") + size(references.ref_fasta, "GB") + 5)

  command <<<
    R --vanilla --no-save <<RSCRIPT

    library(vcfR)
    library(simuscopR)

    if("~{library_type}" == "exome"){

      system(paste0("bamtobed -i ~{emp_bam} > bed_file"))

      seqToProfile("~{emp_bam}", "bed_file", "~{vcf}",
             "~{references.ref_fasta}", "profile")

    } else {
      seqToProfile("~{emp_bam}", vcf.file =  "~{vcf}",
             reference = "~{references.ref_fasta}",  out.profile = "sample.profile")
    }

    RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/reads2map:0.0.1"
    memory: "3 GB"
    cpu: 1
    preemptible: 3
    disks: "local-disk " + disk_size + " HDD"
  }

  output{
    File profile = "sample.profile"
  }
}

task SimuscopSimulation{
 input {
    String library_type
    String sampleName
    Int depth
    File? emp_bam
    File vcf
    Reference references
    String chrom
    File profile
  }

  Int disk_size = ceil(size(vcf, "GB") + size(references.ref_fasta, "GB") + 5)

  command <<<
    R --vanilla --no-save <<RSCRIPT
      vcfR.object <- read.vcfR("~{vcf}")

      variants <- vcf2variants(vcfR.object, sample = "~{sampleName}", chrom = "~{chrom}")

      write.table(variants$SNVs, file = "SNVs.txt", sep = "\t", quote = F, col.names = F, row.names = F)
      write.table(variants$indels, file = "indels.txt", sep = "\t", quote = F, col.names = F, row.names = F)
      write.table(variants$insertions, file = "insertions.txt", sep = "\t", quote = F, col.names = F, row.names = F)

      system("cat SNVs.txt indels.txt insertions.txt > variants.txt")

      if("~{library_type}" == "exome"){
        simuReads(ref = "~{references.ref_fasta}",
              profile = "~{profile}",
              variation = "variants.txt",
              target = "bed_file",
              name = "~{sampleName}",
              output = ".",
              layout = "SE",
              threads = 6,
              verbose = 1,
              coverage = ~{depth})
      } else {
        simuReads(ref = "~{references.ref_fasta}",
              profile = "profile",
              variation = "variants.txt",
              name = "~{sampleName}",
              output = ".",
              layout = "SE", # only single-end by now
              threads = 6,
              verbose = 1,
              coverage = ~{depth})
      }

    RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/reads2map:0.0.1"
    memory: "8 GB"
    cpu:1
    preemptible: 3
    disks: "local-disk " + disk_size + " HDD"
   }

  output {
    File fastq_seq = "~{sampleName}.fq"
  }
}


task RADinitioSimulation{
  input {
    File simu_vcf
    File radinitio_vcf
    String enzyme1
    String? enzyme2
    Reference references
    Int depth
    Int depth_parents
    Int? insert_size
    Int? insert_size_dev
    Int? pcr_cycles
    Int? read_length
    String library_type
    Array[String] names
    String chrom
  }

  # Difficult to guess how much disk we'll need here.
  Int disk_size = ceil(size(simu_vcf, "GB") + size(radinitio_vcf, "GB") + size(references.ref_fasta, "GB") + 5)

  command <<<

    echo -e ~{sep=" " names} > temp
    echo "~{chrom}" > chrom.list

    tr -s ' ' '\n' < temp > temp2
    sed 's/$/\tpop0/' temp2 > popmap.tsv

    mkdir simu_inputs_progeny simu_inputs_parents \
          results_progeny results_parents

    mkdir simu_inputs_progeny/msprime_vcfs simu_inputs_progeny/ref_loci_vars \
          simu_inputs_parents/msprime_vcfs simu_inputs_parents/ref_loci_vars

    # Separate progeny from parents because of the different depths
    head -n 2 popmap.tsv > simu_inputs_parents/popmap.tsv
    lines=$(wc -l popmap.tsv | cut -f1 -d' ')
    tail -n $((lines -2)) popmap.tsv > simu_inputs_progeny/popmap.tsv

    vcftools --vcf ~{simu_vcf} --indv P1 --indv P2 --recode --out parents
    vcftools --vcf ~{simu_vcf} --remove-indv P1 --remove-indv P2 --recode --out progeny

    vcftools --vcf ~{radinitio_vcf} --indv P1 --indv P2 --recode --out parents.rad
    vcftools --vcf ~{radinitio_vcf} --remove-indv P1 --remove-indv P2 --recode --out progeny.rad

    gzip parents.recode.vcf
    gzip parents.rad.recode.vcf
    mv parents.rad.recode.vcf.gz simu_inputs_parents/msprime_vcfs/~{chrom}.vcf.gz
    mv parents.recode.vcf.gz simu_inputs_parents/ref_loci_vars/ri_master.vcf.gz

    gzip progeny.recode.vcf
    gzip progeny.rad.recode.vcf
    mv progeny.rad.recode.vcf.gz simu_inputs_progeny/msprime_vcfs/~{chrom}.vcf.gz
    mv progeny.recode.vcf.gz simu_inputs_progeny/ref_loci_vars/ri_master.vcf.gz


    # progeny
    radinitio --make-library-seq \
              --genome ~{references.ref_fasta} \
              --chromosomes chrom.list \
              --out-dir results_progeny/ \
              --make-pop-sim-dir simu_inputs_progeny/ \
              --library-type ~{library_type} \
              --enz ~{enzyme1} \
              ~{"--enz2 " + enzyme2} \
              --insert-mean ~{default="350" insert_size} \
              --insert-stdev ~{default="35" insert_size_dev} \
              --pcr-cycles ~{default="9" pcr_cycles} \
              --coverage ~{default="20" depth} \
              --read-length ~{default="150" read_length}

    # parents
    radinitio --make-library-seq \
          --genome ~{references.ref_fasta} \
          --chromosomes chrom.list \
          --out-dir results_parents/ \
          --make-pop-sim-dir simu_inputs_parents/ \
          --library-type ~{library_type} \
          --enz ~{enzyme1} \
          ~{"--enz2 " + enzyme2} \
          --insert-mean ~{default="350" insert_size} \
          --insert-stdev ~{default="35" insert_size_dev} \
          --pcr-cycles ~{default="9" pcr_cycles} \
          --coverage ~{default="20" depth_parents} \
          --read-length ~{default="150" read_length}

    # Add fake phred score of 40 (H in Illumina 1.8+ Phred+33)
    # Only in forward read
    for i in results_progeny/rad_reads/*.1.fa.gz; do /seqtk/./seqtk seq -F 'I' $i > $(basename ${i/.fa.gz}.fq); done
    for i in results_parents/rad_reads/*.1.fa.gz; do /seqtk/./seqtk  seq -F 'I' $i > $(basename ${i/.fa.gz}.fq); done

  >>>

  runtime{
    docker: "cristaniguti/radinitio:0.0.1"
    memory: "3 GB"
    cpu:1
    preemptible: 3
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    Array[File] fastq_rad = glob("*.fq")
  }
}
