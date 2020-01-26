version 1.0

import "../structs/reads_simuS.wdl"

workflow reads_simu{

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

  call PedigreeSim2vcf {
    input:
      genotypes_dat = RunPedigreeSimulator.genotypes_dat,
      map_file      = CreatePedigreeSimulatorInputs.mapfile,
      chrom_file    = CreatePedigreeSimulatorInputs.chromfile,
      tot_mks       = CreatePedigreeSimulatorInputs.tot_mks
  }

  call GenerateSampleNames {
    input:
      simulated_vcf = PedigreeSim2vcf.simu_vcf
  }

  scatter (sampleName in GenerateSampleNames.names) {

    call RunVcf2diploid {
      input:
        sampleName = sampleName,
        ref_genome = references.ref_fasta,
        simu_vcf   = PedigreeSim2vcf.simu_vcf
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

    call RunBwaAlignment {
      input:
        sampleName = sampleName,
        reads1     = SimulateIlluminaReads.reads1,
        ref        = references.ref_fasta,
        geno_amb   = references.ref_amb,
        geno_ann   = references.ref_ann,
        geno_bwt   = references.ref_bwt,
        geno_pac   = references.ref_pac,
        geno_sa    = references.ref_sa
    }

    call AddAlignmentHeader {
      input:
        sampleName = sampleName,
        bam_file   = RunBwaAlignment.bam_file,
        bam_idx    = RunBwaAlignment.bam_idx
    }

    call HaplotypeCallerERC {
      input:
        ref        = references.ref_fasta,
        geno_fai   = references.ref_fasta_index,
        sampleName = sampleName,
        bam_rg     = AddAlignmentHeader.bam_rg,
        bam_rg_idx = AddAlignmentHeader.bam_rg_index,
        geno_dict  = references.ref_dict
    }
  }

  call CreateGatkDatabase {
    input:
      path_gatkDatabase = "my_database",
      GVCFs             = HaplotypeCallerERC.GVCF,
      GVCFs_idx         = HaplotypeCallerERC.GVCF_idx
  }

  call GenotypeGVCFs {
    input:
      workspace_tar       = CreateGatkDatabase.workspace_tar,
      output_vcf_filename = "gatk.vcf",
      ref                 = references.ref_fasta,
      geno_fai            = references.ref_fasta_index,
      geno_dict           = references.ref_dict
  }

  call RunFreebayes {
    input:
      freebayesVCFname = "freebayes.vcf",
      ref              = references.ref_fasta,
      bam_rg           = AddAlignmentHeader.bam_rg
  }

  call VcftoolsApplyFilters{
    input:
      freebayesVCF = RunFreebayes.freebayesVCF,
      gatkVCF      = GenotypeGVCFs.gatkVCF
  }

  call CalculateVcfMetrics {
    input:
      freebayesVCF  = VcftoolsApplyFilters.freebayesVCF_F,
      gatkVCF       = VcftoolsApplyFilters.gatkVCF_F,
      tot_mks       = CreatePedigreeSimulatorInputs.tot_mks,
      maternal_trim = SimulateRADseq.maternal_trim,
      seed          = family.seed,
      depth         = family.depth
  }

  Array[Pair[File, String]] bams_files = zip(AddAlignmentHeader.bam_rg, GenerateSampleNames.names)

  scatter (bams in bams_files) {

    call BamCounts{
      input:
        sampleName     = bams.right,
        bam_file       = bams.left,
        bam_idx        = AddAlignmentHeader.bam_rg_index,
        ref            = references.ref_fasta,
        ref_fai        = references.ref_fasta_index,
        ref_dict       = references.ref_dict,
        gatk_vcf       = VcftoolsApplyFilters.gatkVCF_F,
        freebayes_vcf  = VcftoolsApplyFilters.freebayesVCF_F
    }
  }

  call BamCounts4Onemap{
    input:
      sampleName       = GenerateSampleNames.names,
      freebayes_counts = BamCounts.freebayes_counts,
      gatk_counts      = BamCounts.gatk_counts,
      freebayes_pos    = CalculateVcfMetrics.freebayes_pos,
      gatk_pos         = CalculateVcfMetrics.gatk_pos
  }

  Array[String] methods                     = ["gatk", "freebayes"]
  Array[File] vcfs                          = [VcftoolsApplyFilters.gatkVCF_F, VcftoolsApplyFilters.freebayesVCF_F]
  Array[Pair[String, File]] program_and_vcf = zip(methods, vcfs)

  scatter (vcf in program_and_vcf){
    call CreateMaps{
      input:
        tot_mks                   = CreatePedigreeSimulatorInputs.tot_mks,
        simu_vcf                  = PedigreeSim2vcf.simu_vcf,
        methodName                = vcf.left,
        vcf_file                  = vcf.right,
        freebayes_ref_depth       = BamCounts4Onemap.freebayes_ref_bam,
        freebayes_alt_depth       = BamCounts4Onemap.freebayes_alt_bam,
        gatk_ref_depth            = BamCounts4Onemap.gatk_ref_bam,
        gatk_alt_depth            = BamCounts4Onemap.gatk_alt_bam,
        gatk_example_alleles      = BamCounts4Onemap.gatk_example_alleles,
        freebayes_example_alleles = BamCounts4Onemap.freebayes_example_alleles,
        cross                     = family.cross,
        real_phases               = CreatePedigreeSimulatorInputs.real_phases,
        cmBymb                    = family.cmBymb
    }
  }

  call CreateTables{
    input:
        depth                     = family.depth,
        seed                      = family.seed,
        tot_mks                   = CreatePedigreeSimulatorInputs.tot_mks,
        gatk_ref_depth            = CalculateVcfMetrics.gatk_ref_depth,
        gatk_ref_depth_bam        = BamCounts4Onemap.gatk_ref_bam,
        gatk_alt_depth            = CalculateVcfMetrics.gatk_alt_depth,
        gatk_alt_depth_bam        = BamCounts4Onemap.gatk_alt_bam,
        freebayes_ref_depth_bam   = BamCounts4Onemap.freebayes_ref_bam,
        freebayes_alt_depth_bam   = BamCounts4Onemap.freebayes_alt_bam,
        freebayes_ref_depth       = CalculateVcfMetrics.freebayes_ref_depth,
        freebayes_alt_depth       = CalculateVcfMetrics.freebayes_alt_depth,
        all_maps                  = CreateMaps.all_maps,
        all_errors                = CreateMaps.all_errors,
        all_filters               = CreateMaps.all_filters,
        times                     = CreateMaps.times,
        all_RDatas                = CreateMaps.all_RDatas
    }

  output {
    File data1_depths_geno_prob   = CreateTables.data1_depths_geno_prob
    File data2_maps               = CreateTables.data2_maps
    File data3_filters            = CreateTables.data3_filters
    File data5_SNPcall_efficiency = CalculateVcfMetrics.data5_SNPcall_efficiency
    File data4_times              = CreateTables.data4_times
    File data6_RDatas             = CreateTables.data6_RDatas
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
    sed -i 's+chromosome.txt+~{chromfile}+g' ~{parfile}
    sed -i 's+mapfile.txt+~{mapfile}+g' ~{parfile}
    sed -i 's+founderfile.txt+~{founderfile}+g' ~{parfile}
    java -jar /usr/jars/PedigreeSim.jar ~{parfile}
  >>>

  runtime {
    docker: "taniguti/java-in-the-cloud"
  }

  output {
    File genotypes_dat = "sim_genotypes.dat"
  }

}

# Parse pedsim output (.dat) into VCF
task PedigreeSim2vcf {

  input {
    File genotypes_dat
    File map_file
    File chrom_file
    File tot_mks
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT

    library(onemap)
    mks <- read.table("~{tot_mks}", stringsAsFactors = FALSE)
    pos <- mks[,2]
    chr <- mks[,1]

    pedsim2vcf(inputfile = "~{genotypes_dat}",
      map.file = "~{map_file}",
      chrom.file = "~{chrom_file}",
      out.file = "simu.vcf",
      miss.perc = 0, counts = FALSE,pos = pos, haplo.ref = "P1_1",
      chr = chr, phase = TRUE)

    RSCRIPT
  >>>

  runtime {
    docker: "taniguti/onemap"
  }

  output {
    File simu_vcf = "simu.vcf"
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
      --threads=2 \
      --base-calling-profile=~{base_calling} \
      --indel-error-profile=~{indel_error} \
      --gc-bias-profile=~{gc_bias}

  >>>

  runtime {
    docker: "taniguti/pirs-ddrad-cutadapt"
    maxRetries: 3
  }

  output {
    File reads1 = "${sampleName}_100_150_1.fq.gz"
    File reads2 = "${sampleName}_100_150_2.fq.gz"
  }
}

task RunBwaAlignment {

  input {
    String sampleName
    File ref
    File reads1
    File geno_amb
    File geno_ann
    File geno_bwt
    File geno_pac
    File geno_sa
  }

  command <<<
    export PATH=$PATH:/bin
    export PATH=$PATH:/picard.jar
    bwa mem ~{ref} ~{reads1}  | \
    java -jar /picard.jar SortSam \
    I=/dev/stdin \
    O=~{sampleName}.sorted.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true
    mv ~{sampleName}.sorted.bai ~{sampleName}.sorted.bam.bai
  >>>

  runtime {
    docker: "kfdrc/bwa-picard:latest-dev"
  }

  output {
    File bam_file = "${sampleName}.sorted.bam"
    File bam_idx = "${sampleName}.sorted.bam.bai"
  }
}

# Add info to alignment header
task AddAlignmentHeader {
  input {
    String sampleName
    File bam_file
    File bam_idx
  }

  command <<<
    mkdir tmp
    java -jar /gatk/picard.jar AddOrReplaceReadGroups \
      I=~{bam_file} \
      O=~{sampleName}_rg.bam \
      RGLB=lib-~{sampleName} \
      RGPL=illumina \
      RGID=FLOWCELL1.LANE1.~{sampleName} \
      RGSM=~{sampleName} \
      RGPU=FLOWCELL1.LANE1.~{sampleName} \
      CREATE_INDEX=true \
      TMP_DIR=tmp

    mv ~{sampleName}_rg.bai ~{sampleName}_rg.bam.bai
  >>>

  runtime {
    docker: "taniguti/gatk-picard"
  }

  output {
    File bam_rg = "${sampleName}_rg.bam"
    File bam_rg_index = "${sampleName}_rg.bam.bai"
  }
}

# GATK to generate gVCF with variants
task HaplotypeCallerERC {
  input {
    File ref
    File geno_fai
    String sampleName
    File bam_rg
    File bam_rg_idx
    File geno_dict
  }

  command <<<
    /gatk/gatk HaplotypeCaller \
      -ERC GVCF \
      -R ~{ref} \
      -I ~{bam_rg} \
      -O ~{sampleName}_rawLikelihoods.g.vcf \
      --max-reads-per-alignment-start 0
  >>>

  runtime {
    docker: "taniguti/gatk-picard"
  }

  output {
    File GVCF = "${sampleName}_rawLikelihoods.g.vcf"
    File GVCF_idx = "${sampleName}_rawLikelihoods.g.vcf.idx"
  }
}

task CreateGatkDatabase {

  input {
    String path_gatkDatabase
    Array[File] GVCFs
    Array[File] GVCFs_idx
  }

  command <<<
    /gatk/gatk GenomicsDBImport \
      --genomicsdb-workspace-path ~{path_gatkDatabase} \
      -L Chr10 \
      -V ~{sep=" -V "  GVCFs}

    tar -cf ~{path_gatkDatabase}.tar ~{path_gatkDatabase}
  >>>

  runtime {
    docker: "taniguti/gatk-picard"
  }

  output {
    File workspace_tar = "${path_gatkDatabase}.tar"
  }
}

# Variant calling on gVCF
task GenotypeGVCFs {

  input {
    File workspace_tar
    String output_vcf_filename
    File ref
    File geno_fai
    File geno_dict
  }

  command <<<
    tar -xf ~{workspace_tar}
    WORKSPACE=$( basename ~{workspace_tar} .tar)

    /gatk/gatk GenotypeGVCFs \
        -R ~{ref} \
        -O ~{output_vcf_filename} \
        -G StandardAnnotation \
        -V gendb://$WORKSPACE
  >>>

  runtime {
    docker: "taniguti/gatk-picard"
  }

  output {
    File gatkVCF = "${output_vcf_filename}"
    File gatkVCF_index = "${output_vcf_filename}.idx"
  }

}

# Variant calling using freebayes
task RunFreebayes {

  input {
    String freebayesVCFname
    File ref
    Array[File] bam_rg
  }

  command <<<
    freebayes --genotype-qualities -f ~{ref} ~{sep=" "  bam_rg} > ~{freebayesVCFname}
  >>>

  runtime {
    docker: "taniguti/freebayes"
  }

  output {
    File freebayesVCF = "${freebayesVCFname}"
  }
}


task  VcftoolsApplyFilters {

  input {
    File gatkVCF
    File freebayesVCF
  }

  command <<<
    vcftools --vcf "~{gatkVCF}" --max-missing 0.75 --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --out gatk
    vcftools --vcf "~{freebayesVCF}" --max-missing 0.75 --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --out freebayes

  >>>
  runtime {
    docker: "taniguti/vcftools"
  }

  output {
    File gatkVCF_F = "gatk.recode.vcf"
    File freebayesVCF_F = "freebayes.recode.vcf"
  }
}

task CalculateVcfMetrics {

  input {
    File freebayesVCF
    File gatkVCF
    File tot_mks
    Array[File] maternal_trim
    Int seed
    Int depth
  }

  command <<<

        R --vanilla --no-save <<RSCRIPT

        system("cp ~{sep=" "  maternal_trim} .")
        library(vcfR)
        freebayes <- read.vcfR("~{freebayesVCF}")
        gatk <- read.vcfR("~{gatkVCF}")
        maternal <- strsplit("~{sep=" ; "  maternal_trim}", split=";")[[1]][1]
        system(paste("grep '>'", maternal ,"> frags"))

        frags <- read.table("frags", stringsAsFactors=F)
        start <- frags[,14]
        end <- start + 202

        snps <- read.table("~{tot_mks}", stringsAsFactors = F)
        real.pos <- snps[,2]

        filt.idx <- vector()
        for(i in 1:length(start))
        filt.idx <- c(filt.idx,which(real.pos >= start[i] & real.pos <= end[i]))

        snps.filt <- snps[filt.idx,]
        filt.pos <- snps.filt[,2]
        ref.filt <- snps.filt[,3]
        alt.filt <- snps.filt[,4]

        methods <- c("freebayes", "gatk") # include in a scatter
        results_tot <- vector()
        for(i in methods){

          # counting corrected identified markers
          pos <- get(i)@fix[,2]
          chr <- get(i)@fix[,1]
          site_list <- data.frame(chr, pos, pos)
          # Export for next step
          write.table(site_list, file= paste0(i,"_site_list.txt"), quote=F, row.names=F, sep="\t", col.names=F)

          ref <- get(i)@fix[,4]
          alt <- get(i)@fix[,5]

          nmk.filt <- length(filt.pos)
          nmk.id <- length(pos)

          ok <- sum(filt.pos %in% pos) #  marcadores identificados do total
          falso.positivo <- sum(!(pos %in% filt.pos)) # falsos positivos
          ref.ok <- sum(ref.filt==ref[pos %in% filt.pos])
          alt.ok <- sum(alt.filt==alt[pos %in% filt.pos])

          result <- data.frame(depth = ~{depth}, seed = ~{seed}, SNPcall = i,mks_tot = nmk.filt, mks_ide = nmk.id, ok, fake=falso.positivo, ref.ok, alt.ok)
          results_tot <- rbind(results_tot, result)

          #write.table(result, file= paste0(i,".txt"), quote=F, row.names=F, sep="\t")

          # tables for mesure depth distribuition
          if(dim(get(i)@gt)[1] != 0){
              idx <- which(strsplit(get(i)@gt[1,1], split=":")[[1]] == "AD")
              if(length(idx) != 0){
              ref.depth <- matrix(sapply(strsplit(sapply(strsplit(get(i)@gt[,-1], split=":"), "[",idx), split=","), "[",1),
                                  ncol = dim(get(i)@gt)[2]-1)
              alt.depth <- matrix(sapply(strsplit(sapply(strsplit(get(i)@gt[,-1], split=":"), "[",idx), split=","), "[",2),
                                  ncol =  dim(get(i)@gt)[2]-1)
              colnames(ref.depth) <- colnames(alt.depth) <- colnames(get(i)@gt[,-1])
              rownames(ref.depth) <- rownames(alt.depth) <- paste0(get(i)@fix[,1],"_", get(i)@fix[,2])
              write.table(ref.depth, file = paste0(i,"_ref_depth.txt"), quote=F, row.names=T, sep="\t")
              write.table(alt.depth, file = paste0(i, "_alt_depth.txt"), quote=F, row.names=T, sep="\t")
              } else {
                  null.table <- matrix(rep(0,dim(get(i)@gt)[2]-1), ncol=dim(get(i)@gt)[2]-1)
                  write.table(null.table, file = paste0(i,"_ref_depth.txt"), quote=F, row.names=F, sep="\t")
                  write.table(null.table, file = paste0(i, "_alt_depth.txt"), quote=F, row.names=F, sep="\t")
              }
              # table for GQ
              idx <- which(strsplit(get(i)@gt[1,1], split=":")["FORMAT"] == "GQ")
              if(length(idx)!=0){
                  GQ <- sapply(strsplit(get(i)@gt[,-1], split=":"), "[",idx)
                  write.table(GQ, file = paste0(i, "_GQ.txt"), quote=F, row.names=F, sep="\t")
                  } else {
                      null.table <- matrix(rep(0,dim(get(i)@gt)[2]-1), ncol=dim(get(i)@gt)[2]-1)
                      write.table(null.table, file = paste0(i, "_GQ.txt"), quote=F, row.names=F, sep="\t")
                  }

              } else{
                  null.table <- matrix(rep(0,dim(get(i)@gt)[2]-1), ncol=dim(get(i)@gt)[2]-1)
                  write.table(null.table, file = paste0(i,"_ref_depth.txt"), quote=F, row.names=F, sep="\t")
                  write.table(null.table, file = paste0(i, "_alt_depth.txt"), quote=F, row.names=F, sep="\t")
                  write.table(null.table, file = paste0(i, "_GQ.txt"), quote=F, row.names=F, sep="\t")
              }
        }

        saveRDS(results_tot, file= "data5_SNPcall_efficiency.rds")
        RSCRIPT

  >>>

  runtime {
    docker: "taniguti/onemap"
  }

  output {
    File freebayes_pos = "freebayes_site_list.txt"
    File gatk_pos = "gatk_site_list.txt"
    File data5_SNPcall_efficiency  = "data5_SNPcall_efficiency.rds"
    File freebayes_ref_depth = "freebayes_ref_depth.txt"
    File freebayes_alt_depth = "freebayes_alt_depth.txt"
    File gatk_ref_depth = "gatk_ref_depth.txt"
    File gatk_alt_depth = "gatk_alt_depth.txt"
  }
}

# This task extract the allele depths from bam files
task BamCounts{
  input{
    String sampleName
    File bam_file
    Array[File] bam_idx
    File ref
    File ref_fai
    File ref_dict
    File gatk_vcf
    File freebayes_vcf
  }

  command <<<

    java -jar /gatk/picard.jar VcfToIntervalList \
      I=~{gatk_vcf} \
      O=gatk.interval.list

    /gatk/gatk CollectAllelicCounts \
      --input ~{bam_file} \
      --reference ~{ref} \
      --intervals gatk.interval.list \
      --output ~{sampleName}_gatk_counts.tsv

   java -jar /gatk/picard.jar VcfToIntervalList \
      I=~{freebayes_vcf} \
      O=freebayes.interval.list

    /gatk/gatk CollectAllelicCounts \
      --input ~{bam_file} \
      --reference ~{ref} \
      --intervals freebayes.interval.list \
      --output ~{sampleName}_freebayes_counts.tsv

  >>>

  runtime{
    docker:"taniguti/gatk-picard"
  }

  output{
    File gatk_counts = "~{sampleName}_gatk_counts.tsv"
    File freebayes_counts = "~{sampleName}_freebayes_counts.tsv"
  }
}

# This task convert the output from BamCounts to the depths input for onemap
task BamCounts4Onemap{
  input{
    Array[File] freebayes_counts
    Array[File] gatk_counts
    Array[String] sampleName
    File freebayes_pos
    File gatk_pos
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(R.utils)
      system("cp ~{sep=" "  freebayes_counts} .")
      system("cp ~{sep=" "  gatk_counts} .")
      system("cp ~{freebayes_pos} .")
      system("cp ~{gatk_pos} .")
      names <- c("~{sep=" , "  sampleName}")
      names <- unlist(strsplit(names, split = " , "))

      methods <- c("gatk", "freebayes")

      for(method in methods){

          file.counts <- read.table(paste0(names[1],"_", method,"_counts.tsv"), skip = 3, header=T, stringsAsFactors = F)

          ref_depth_matrix2 <- alt_depth_matrix2  <- matrix(NA, nrow = dim(file.counts)[1], ncol = length(names))

          for(j in 1:length(names)){
            ## From picard tool

            file.counts <- read.table(paste0(names[j],"_", method,"_counts.tsv"), skip = 3, header=T, stringsAsFactors = F)

            ref_depth_matrix2[,j] <- file.counts[,3]
            alt_depth_matrix2[,j] <- file.counts[,4]

            if(j == 1){
              ref_allele <- file.counts[,5]
              alt_allele <- file.counts[,6]
            } else {
              idx.ref <- which(ref_allele == "N")
              idx.alt <- which(alt_allele == "N")
              if(length(idx.ref)!=0){
                ref_allele[idx.ref] <- file.counts[idx.ref,5]
              }
              if(length(idx.alt)!=0){
                alt_allele[idx.alt] <- file.counts[idx.alt,6]
              }
            }

          }

          rownames(ref_depth_matrix2) <- rownames(alt_depth_matrix2) <- paste0(file.counts[,1],"_", file.counts[,2])
          colnames(ref_depth_matrix2) <- colnames(alt_depth_matrix2) <- names

          alleles <- data.frame(file.counts[,1],file.counts[,2], ref_allele, alt_allele)
          write.table(alleles, file = paste0(method,"_example4ref_alt_alleles.txt"), col.names = F, row.names = F)

          write.table(ref_depth_matrix2, file = paste0(method,"_ref_depth_bam.txt"), quote=F, row.names=T, sep="\t", col.names=T)
          write.table(alt_depth_matrix2, file = paste0(method,"_alt_depth_bam.txt"), quote=F, row.names=T, sep="\t", col.names=T)
        }

    RSCRIPT
  >>>

  runtime{
    docker:"taniguti/onemap"
  }

  output{
    File gatk_ref_bam      = "gatk_ref_depth_bam.txt"
    File gatk_alt_bam      = "gatk_alt_depth_bam.txt"
    File freebayes_ref_bam = "freebayes_ref_depth_bam.txt"
    File freebayes_alt_bam = "freebayes_alt_depth_bam.txt"
    File gatk_example_alleles    = "gatk_example4ref_alt_alleles.txt"
    File freebayes_example_alleles    = "freebayes_example4ref_alt_alleles.txt"

  }
}

task CreateMaps {

 input {
    File tot_mks
    String methodName
    File simu_vcf
    File vcf_file
    File freebayes_ref_depth
    File freebayes_alt_depth
    File gatk_ref_depth
    File gatk_alt_depth
    File gatk_example_alleles
    File freebayes_example_alleles
    String cross
    File real_phases
    Float cmBymb
  }

  command <<<

        R --vanilla --no-save <<RSCRIPT
          library(supermassa4onemap)
          library(onemap)
          library(updog)
          library(reshape2)
          library(vcfR)
          library(doParallel)
          library(GUSMap)

          args = commandArgs(trailingOnly=TRUE)

          system("cp ~{freebayes_ref_depth} .")
          system("cp ~{freebayes_alt_depth} .")
          system("cp ~{gatk_ref_depth} .")
          system("cp ~{gatk_alt_depth} .")
          system("cp ~{gatk_example_alleles} .")
          system("cp ~{freebayes_example_alleles} .")
          method_name <- "~{methodName}"
          tot_mks_file <- "~{tot_mks}"
          simu_vcf_file <- "~{simu_vcf}"
          vcf_file <- "~{vcf_file}"
          cross <- "~{cross}"
                cMbyMb <- ~{cmBymb}
          real_phases <- read.table("~{real_phases}")
          source("/opt/scripts/functions_simu.R")


          ## KNOWN VARIANTS
          tot_mks <- read.table(tot_mks_file)
          
          if(cross == "F1"){
            cross <- "outcross"
            f1 = NULL
          } else if (cross == "F2"){
            cross <- "f2 intercross"
            f1 = "F1"
          }
          
          # READING DATA FROM SIMULATED POPULATION
          simu <- read.vcfR(simu_vcf_file)
          gab <- onemap_read_vcfR(vcfR.object=simu,
                                  cross= cross,
                                  parent1="P1",
                                  parent2="P2",
                                  f1 = f1)
          
          ## READING FINAL VCF FROM PIPELINE
          vcf <- read.vcfR(vcf_file)
          df <- onemap_read_vcfR(vcfR.object=vcf,
                                 cross= cross,
                                 parent1="P1",
                                 parent2="P2",
                                 f1 = f1)
          
          ## FILTERS REPORT
          SNPcall <- method_name
          Genocall <- "df"
          CountsFrom <- "vcf"
          filters_tab <- create_filters_report(df, SNPcall, CountsFrom, Genocall)
          
          ## MAPS REPORT - DF
          
          times <-system.time(create_maps_report(input.seq = filters_tab, 
                                                 tot_mks = tot_mks, gab = gab, 
                                                 SNPcall , Genocall,
                                                 fake= F, CountsFrom,cMbyMb))
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", FALSE)
          times <- data.frame(meth = outname, time = times[3])
          
          times_temp <-system.time(create_maps_report(input.seq = filters_tab, 
                                                      tot_mks = tot_mks, gab = gab, 
                                                      SNPcall = method_name, Genocall= "df",
                                                      fake= T, CountsFrom="vcf",cMbyMb))
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", TRUE)
          times_temp <- data.frame(meth = outname, time = times_temp[3])
          times <- rbind(times, times_temp)
          
          errors_tab <- create_errors_report(onemap_obj = df, gab,
                                             SNPcall, Genocall,
                                             CountsFrom)
          
          # MAPS REPORT - GQ
          aval.gq <- extract_depth(vcfR.object=vcf,
                                   onemap.object=df,
                                   vcf.par="GQ",
                                   parent1="P1",
                                   parent2="P2",
                                   f1 = f1,
                                   recovering=FALSE)
          
          aval.gq <- create_probs(df, genotypes_errors=aval.gq)
          
          Genocall <- "GQ"
          filters_tab <- create_filters_report(aval.gq, SNPcall, CountsFrom, Genocall)
          
          fake <- F
          
          times_temp <- system.time(create_maps_report(filters_tab, 
                                                       tot_mks, gab, 
                                                       SNPcall, Genocall,
                                                       fake, CountsFrom,cMbyMb))
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
          times_temp <- data.frame(meth = outname, time = times_temp[3])
          times <- rbind(times, times_temp)
          
          fake <- T
          times_temp <- system.time(create_maps_report(filters_tab, tot_mks, gab, 
                                                       SNPcall, Genocall,
                                                       fake, CountsFrom,cMbyMb))
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
          times_temp <- data.frame(meth = outname, time = times_temp[3])
          times <- rbind(times, times_temp)
          
          
          create_errors_report(aval.gq, gab,
                               SNPcall, Genocall,
                               CountsFrom)
          
          
          # OTHER TOOLS
          ## With depths from vcf
          
          updog.aval <- updog_error(
            vcfR.object=vcf,
            onemap.object=df,
            vcf.par="AD",
            parent1="P1",
            parent2="P2",
            f1 = f1,
            recovering=TRUE,
            mean_phred=20,
            cores=3,
            depths=NULL)
          
          supermassa.aval <- supermassa4onemap::supermassa_error(
            vcfR.object=vcf,
            onemap.object = df,
            vcf.par = "AD",
            parent1 = "P1",
            parent2 = "P2",
            f1 = f1,
            recovering = TRUE,
            mean_phred = 20,
            cores = 3,
            depths = NULL)
          
          polyrad.aval <- polyRAD_error(
            vcf=vcf_file,
            onemap.obj=df,
            parent1="P1",
            parent2="P2",
            f1 = f1,
            crosstype=cross)
          
          metodologies <- list(updog = updog.aval, supermassa = supermassa.aval, polyrad = polyrad.aval)
          for (metodology in names(metodologies)){
            error_aval <- metodologies[[metodology]]
            ## Filters
            Genocall <- metodology
            SNPcall <- method_name
            CountsFrom <- "vcf"
            
            filters_tab <- create_filters_report(error_aval, SNPcall, CountsFrom, Genocall)
            
            ## Maps
            fake <- F
            times_temp <- system.time(create_maps_report(input.seq = filters_tab, 
                                                         tot_mks, gab, 
                                                         SNPcall, Genocall,
                                                         fake, CountsFrom,cMbyMb))
            
            outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
            times_temp <- data.frame(meth = outname, time = times_temp[3])
            times <- rbind(times, times_temp)
            
            fake <- T
            times_temp <- system.time(create_maps_report(input.seq = filters_tab, 
                                                         tot_mks, gab, 
                                                         SNPcall, Genocall,
                                                         fake, CountsFrom,cMbyMb))
            
            outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
            times_temp <- data.frame(meth = outname, time = times_temp[3])
            times <- rbind(times, times_temp)
            
            ## Errors info
            
            create_errors_report(error_aval, gab,  
                                 SNPcall, Genocall,
                                 CountsFrom)
            
          }
          
          
          ## Depths from bam
          depths.alt <- read.table(paste0(method_name, "_alt_depth_bam.txt"), header = T)
          depths.ref <- read.table(paste0(method_name, "_ref_depth_bam.txt"), header = T)
          
          depths <- list("ref"=depths.ref, "alt"=depths.alt)
          CountsFrom <- "bam"
          updog.aval.bam <- updog_error(
            vcfR.object=vcf,
            onemap.object=df,
            vcf.par="AD",
            parent1="P1",
            parent2="P2",
            f1 = f1,
            recovering=TRUE,
            mean_phred=20,
            cores=3,
            depths=depths)
          
          supermassa.aval.bam <- supermassa_error(
            vcfR.object=vcf,
            onemap.object = df,
            vcf.par = "AD",
            parent1 = "P1",
            parent2 = "P2",
            f1 = f1,
            recovering = TRUE,
            mean_phred = 20,
            cores = 3,
            depths = depths)
          
          new.vcf <- make_vcf(vcf_file, depths, method_name)
          
          polyrad.aval.bam <- polyRAD_error(
            vcf=new.vcf,
            onemap.obj=df,
            parent1="P1",
            parent2="P2",
            f1 = f1,
            crosstype=cross)
          
          metodologies <- list(updog = updog.aval.bam, supermassa= supermassa.aval.bam, polyrad=polyrad.aval.bam)
          for (metodology in names(metodologies)){
            error_aval <- metodologies[[metodology]]
            ## Filters
            Genocall <- metodology
            CountsFrom <- "bam"
            
            filters_tab <- create_filters_report(error_aval, SNPcall, CountsFrom, Genocall)
            
            ## Maps
            fake <- F
            times_temp <- system.time(create_maps_report(filters_tab,  
                                                         tot_mks, gab, 
                                                         SNPcall, Genocall,
                                                         fake, CountsFrom,cMbyMb))
            
            outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
            times_temp <- data.frame(meth = outname, time = times_temp[3])
            times <- rbind(times, times_temp)
            
            fake = T
            times_temp <- system.time(create_maps_report(filters_tab,  
                                                         tot_mks, gab, 
                                                         SNPcall, Genocall,
                                                         fake, CountsFrom,cMbyMb))
            
            outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
            times_temp <- data.frame(meth = outname, time = times_temp[3])
            times <- rbind(times, times_temp)
            
            ## Errors info
            errors_tab <- create_errors_report(onemap_obj = error_aval, gab = gab, 
                                               SNPcall, Genocall,
                                               CountsFrom)
          }
          
          ## Gusmap maps
          Genocall <- "gusmap"
          
          fake <- F
          CountsFrom <- "vcf"
          times_temp <- system.time(create_gusmap_report(vcf_file, gab,SNPcall, Genocall,
                                                         fake, CountsFrom))
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
          times_temp <- data.frame(meth = outname, time = times_temp[3])
          times <- rbind(times, times_temp)
          
          fake <- T
          times_temp <- system.time(create_gusmap_report(vcf_file, gab,SNPcall, Genocall,
                                                         fake, CountsFrom))
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
          times_temp <- data.frame(meth = outname, time = times_temp[3])
          times <- rbind(times, times_temp)
          
          
          CountsFrom <- "bam"
          fake <- F
          times_temp <- system.time(create_gusmap_report(new.vcf, gab,SNPcall, Genocall,
                                                         fake, CountsFrom))
          
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
          times_temp <- data.frame(meth = outname, time = times_temp[3])
          times <- rbind(times, times_temp)
          
          fake <- T
          times_temp <- system.time(create_gusmap_report(new.vcf, gab,SNPcall, Genocall,
                                                         fake, CountsFrom))
          
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
          times_temp <- data.frame(meth = outname, time = times_temp[3])
          times <- rbind(times, times_temp)
          
          write.table(times, paste0(method_name,"_times.txt"))
          
          Genocall <- c("df", "GQ", "updog", "supermassa", "polyrad", "gusmap")
          fake <- c(TRUE, FALSE)
          CountsFrom <- c("vcf", "bam")
          
          all_maps <- all_errors <- all_filters <- data.frame()
          all_RDatas <- list()
          z <- 1
          names_RDatas <- vector()
          for(i in 1:length(Genocall)){
            for(j in 1:length(CountsFrom)){
              for(w in 1:length(fake)){
                if(CountsFrom[j] == "bam" & (Genocall[i] == "df" | Genocall[i] == "GQ")){
                } else {
                  cat(paste0("map_",method_name,"_",CountsFrom[j], "_", Genocall[i], "_",fake[w], ".txt"), "\n")
                  temp_map <- read.table(paste0("map_",method_name,"_",CountsFrom[j], "_", Genocall[i], "_",fake[w], ".txt"), header=T)
                  load(paste0("map_",method_name,"_",CountsFrom[j], "_", Genocall[i], "_",fake[w], ".RData"))
                  all_RDatas[[z]] <- map_df
                  names_RDatas <- c(names_RDatas, paste0("map_",method_name,"_",CountsFrom[j], "_", Genocall[i], "_",fake[w]))
                  z <- z+1
                  all_maps <- rbind(all_maps, temp_map)
                  if(Genocall[i] == "gusmap"){
                
                  } else {
                    cat(paste0("errors_",method_name,"_",CountsFrom[j], "_", Genocall[i], ".txt"), "\n")
                    temp_error <- read.table(paste0("errors_",method_name,"_",CountsFrom[j], "_", Genocall[i],".txt"), header=T)
                    all_errors <- rbind(all_errors, temp_error)
                    cat(paste0("filters_",method_name,"_",CountsFrom[j], "_", Genocall[i],".txt"), "\n")
                    temp_filters <- read.table(paste0("filters_",method_name,"_",CountsFrom[j], "_", Genocall[i],".txt"), header=T)
                    all_filters <- rbind(all_filters, temp_filters)
                  }
                }
              }
            }
          }

          names(all_RDatas) <- names_RDatas
          write.table(all_maps, paste0(method_name,"_maps.txt"))
          write.table(all_errors, paste0(method_name,"_errors.txt"))
          write.table(all_filters, paste0(method_name,"_filters.txt"))
          save(all_RDatas, file= paste0(method_name, "_RDatas.RData"))

        RSCRIPT
  >>>

  runtime{
    docker:"taniguti/onemap"
  }

  output{
    File times = "~{methodName}_times.txt"
    File all_maps = "~{methodName}_maps.txt"
    File all_errors = "~{methodName}_errors.txt"
    File all_filters = "~{methodName}_filters.txt"
    File all_RDatas = "~{methodName}_RDatas.RData"
  }
}

task CreateTables{
  input{
    Int depth
    Int seed
    File tot_mks
    File gatk_ref_depth
    File gatk_ref_depth_bam
    File gatk_alt_depth
    File gatk_alt_depth_bam
    File freebayes_ref_depth_bam
    File freebayes_alt_depth_bam
    File freebayes_ref_depth
    File freebayes_alt_depth
    Array[File] all_maps
    Array[File] all_errors
    Array[File] all_filters
    Array[File] times
    Array[File] all_RDatas
  }

  command <<<

        R --vanilla --no-save <<RSCRIPT

          library(reshape2)
          library(vcfR)
          system("cp ~{sep= " " times } .")
          system("cp ~{sep= " " all_maps } .")
          system("cp ~{sep= " " all_errors } .")
          system("cp ~{sep= " " all_filters } .")
          system("cp ~{sep= " " all_RDatas } .")
         
          tot_mks <- read.table("~{tot_mks}")
          method <- c("gatk", "freebayes")
          meth.bam <-  c("polyrad", "updog", "supermassa")
          meth.filt <- c("dfAndGQ", "polyrad", "updog", "supermassa")
          meth.geno <- c("df", "GQ", "polyrad", "updog", "supermassa")
          seed <- ~{seed}
          depth <- ~{depth}

          # Functions

           joint_depths <- function(depth_matrix_alt, depth_matrix_ref, CountsFrom, SNPcall, depth,seed, genotype_meth){
              depth_matrix <- list(depth_matrix_alt, depth_matrix_ref)
              alle_name <- c("alt", "ref")
              allele <- list()
              for(i in 1:length(depth_matrix)){
                depth2 <- read.table(depth_matrix[[i]], header = T, stringsAsFactors = F)
                depth3 <- as.data.frame(apply(depth2, 2, as.integer))
                depth3 <- cbind(MKS= rownames(depth2), depth3)
                depth3 <- depth3[,sort(colnames(depth3))]
                allele[[i]] <- melt(depth3)
                colnames(allele[[i]]) <- c("mks", "ind", paste0(alle_name[i]))
              }
              alleles <- merge(allele[[1]], allele[[2]])
              alleles <- cbind(seed= seed, depth = depth, "SNPcall" = SNPcall, "CountsFrom" = CountsFrom, alleles)
              alleles[,4] <- as.character(alleles[,4])
              alleles[,5] <- as.character(alleles[,5])
              return(alleles)
            }


            ########################################################################################
            # Table1: GenoCall; mks; ind; SNPcall; CountsFrom; alt; ref; gabGT; methGT; A; AB; BA; B
            ########################################################################################
            df_tot_all <- times <- maps_tot <- filters_tot <- vector()
            tot_RDatas <- list()
            for(i in 1:length(method)){

              if(method[i] == "gatk"){
                alt.depth.bam <- "~{gatk_alt_depth_bam}"
                ref.depth.bam <- "~{gatk_ref_depth_bam}"
                alt.depth <- "~{gatk_alt_depth}"
                ref.depth <- "~{gatk_ref_depth}"
              } else{
                alt.depth.bam <- "~{freebayes_alt_depth_bam}"
                ref.depth.bam <- "~{freebayes_ref_depth_bam}"
                alt.depth <- "~{freebayes_alt_depth}"
                ref.depth <- "~{freebayes_ref_depth}"
              }
  
              ## Depths by bam
              
              df_tot_bam <- joint_depths(depth_matrix_alt = alt.depth.bam, depth_matrix_ref = ref.depth.bam, 
                                         CountsFrom = "bam", SNPcall = method[i], depth = depth, seed = seed, genotype_meth = meth.bam)
              
              ## Depths by softwares
              df_tot <- joint_depths(alt.depth,ref.depth, "vcf", method[i], depth, seed, meth.geno)
              
              df_tot_tot <- rbind(df_tot_bam, df_tot)
              
              chr <- unique(sapply(strsplit(df_tot_tot[,5], "_"), "[",1))
              all_errors <- read.table(paste0(method[i],"_errors.txt"), header = T)
              all_errors[,5] <- paste0(chr, "_",all_errors[,5])
              colnames(all_errors) <- c("SNPcall", "Genocall", "CountsFrom", "ind", "mks", "gabGT", "methGT", "A", "AB", "BA", "B")
              df_meth <- merge(all_errors, df_tot_tot, by = c("SNPcall", "CountsFrom", "ind", "mks"))
              
              df_tot_all <- rbind(df_tot_all, df_meth)
              
              ########################################################
              # Table2: seed; CountsFrom; ErrorProb; SNPcall; MK; rf; phases; real_phases
              ########################################################
            
              map_df <- read.table(paste0(method[i], "_maps.txt"), header = T)
              map_df <- data.frame("seed" = seed, "depth" = depth, map_df)
            
              maps_tot <- rbind(maps_tot, map_df)
              
              ##########################################################################
              # Table3: CountsFrom; seed; SNPcall; GenoCall; n_mks; distorted; redundant
              ##########################################################################
              
              filters_temp <- read.table(paste0(method[i], "_filters.txt"), header = T)
              filters_temp <- data.frame("seed" = seed, "depth" = depth, filters_temp)
              
              filters_tot <- rbind(filters_tot, filters_temp)
              
              ###########################################################################
              # Table4: CountsFrom; seed; SNPcall; GenoCall
              ###########################################################################
              
              CountsFrom <- vector()
              times_temp <- read.table(paste0(method[i], "_times.txt"), stringsAsFactors = F)
              temp <- strsplit(times_temp[,1], "_")
              temp <- do.call(rbind, temp)
              CountsFrom[which(temp[,3] == "bam")] <- "bam"
              CountsFrom[which(temp[,3] != "bam")] <- "vcf"
              SNPcall <- temp[,2]
              temp[which(temp[,3] == "bam"),3] <- temp[which(temp[,3] == "bam"),4]
              Genocall <- temp[,4]
              real.mks <- temp[,5]
              times <- rbind(times, data.frame(CountsFrom, seed, depth, SNPcall, Genocall, real.mks, times = times_temp[,2]))

              ###########################################################################
              # Table6: list of RDatas with name CountsFrom; seed; SNPcall; GenoCall
              ###########################################################################

              load(paste0(method[i], "_RDatas.RData"))
              tot_RDatas <- c(tot_RDatas, all_RDatas)

            }

            names_RDatas <- names(tot_RDatas)
            new_names <- paste0(seed, "_", depth, "_", names_RDatas)
            names(tot_RDatas) <- new_names
            
            saveRDS(df_tot_all, file = "data1_depths_geno_prob.rds")
            saveRDS(maps_tot, file = "data2_maps.rds")
            saveRDS(filters_tot, file = "data3_filters.rds")
            saveRDS(times, file= "data4_times.rds")
            save(tot_RDatas, file = "data6_RDatas.RData")

      RSCRIPT
  >>>

  runtime{
    docker:"taniguti/onemap"
  }

  output{
    File data1_depths_geno_prob = "data1_depths_geno_prob.rds"
    File data2_maps = "data2_maps.rds"
    File data3_filters = "data3_filters.rds"
    File data4_times   = "data4_times.rds"
    File data6_RDatas  = "data6_RDatas.RData"
  }
}

