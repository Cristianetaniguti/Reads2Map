version 1.0

import "./structs/reads_simuS.wdl"

workflow reads_simu{

  input {
    ReferenceFasta references
    Family family
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
        depth         = family.depth
    }

    call RunBwaAlignment {
      input:
        sampleName = sampleName,
        reads1     = SimulateIlluminaReads.reads1,
        reads2     = SimulateIlluminaReads.reads2,
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
      output_vcf_filename = family.name + "_gatk.vcf",
      ref                 = references.ref_fasta,
      geno_fai            = references.ref_fasta_index,
      geno_dict           = references.ref_dict
  }

  call RunFreebayes {
    input:
      freebayesVCFname = family.name + "_freebayes.vcf",
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
      map_file      = CreatePedigreeSimulatorInputs.mapfile,
      maternal_trim = SimulateRADseq.maternal_trim
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

  Array[String] methods = ["gatk", "freebayes"]
  Array[File] vcfs = [VcftoolsApplyFilters.gatkVCF_F, VcftoolsApplyFilters.freebayesVCF_F]
  Array[Pair[String, File]] program_and_vcf = zip(methods, vcfs)

  scatter (vcf in program_and_vcf){
    call all_maps{
      input:
        tot_mks    = CreatePedigreeSimulatorInputs.tot_mks,
        simu_vcf   = PedigreeSim2vcf.simu_vcf,
        methodName = vcf.left,
        vcf_file   = vcf.right,
        freebayes_ref_depth2 = BamCounts4Onemap.freebayes_ref_bam, 
        freebayes_alt_depth2 = BamCounts4Onemap.freebayes_alt_bam,
        gatk_ref_depth2 = BamCounts4Onemap.gatk_ref_bam,
        gatk_alt_depth2 = BamCounts4Onemap.gatk_alt_bam,
        gatk_example_alleles = BamCounts4Onemap.gatk_example_alleles,
        freebayes_example_alleles = BamCounts4Onemap.freebayes_example_alleles,
        cross = family.cross
    }
  }

  output {
    File freebayes_vcf       = VcftoolsApplyFilters.freebayesVCF_F
    File gatk_vcf            = VcftoolsApplyFilters.gatkVCF_F
    File freebayes_aval_vcf  = CalculateVcfMetrics.freebayes_aval_vcf
    File gatk_aval_vcf       = CalculateVcfMetrics.gatk_aval_vcf
    File freebayes_ref_depth = CalculateVcfMetrics.freebayes_ref_depth
    File freebayes_alt_depth = CalculateVcfMetrics.freebayes_alt_depth
    File gatk_ref_depth      = CalculateVcfMetrics.gatk_ref_depth
    File gatk_alt_depth      = CalculateVcfMetrics.gatk_alt_depth
    File tot_mks             = CreatePedigreeSimulatorInputs.tot_mks
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
  
        founder1.df <- matrix(ref.alleles, ncol = ploidy, nrow = length(ref.alleles))
  
        founder2.df <- matrix(ref.alleles, ncol = ploidy, nrow = length(ref.alleles))
  
        idx <- 1:length(ref.alleles)
        for(i in 1:length(doses)){
          size <- round((doses[i]/100)*length(ref.alleles))
          if(i == 1){
            idx.both <- sample(idx, as.numeric(size)*2) # It will not have monomorphic markers
            idx.p1 <- idx.both[1:as.numeric(size)]
            idx.p2 <- idx.both[(as.numeric(size)+1):(as.numeric(size)*2)]
            founder1.df[idx.p1,] <- ref.alleles[idx.p1]
            founder2.df[idx.p2,] <- ref.alleles[idx.p2]
            idx.p1.tot <- idx[-idx.p1]
            idx.p2.tot <- idx[-idx.p2]
          } else if(i == 2){
            idx.p1 <- sample(idx.p1.tot, as.numeric(size)-1)
            idx.p2 <- vector()
            for(w in 1:(as.numeric(size)-1)){
              idx.p2[w] <- sample(idx.p2.tot, 1)
              while(any(idx.p1 %in% idx.p2[w])){
                idx.p2[w] <- sample(idx.p2.tot, 1)
              }
              idx.p2.tot <- idx.p2.tot[-which(idx.p2.tot%in%idx.p2)]
            }
            idx.p1.tot <- idx.p1.tot[-which(idx.p1.tot%in%idx.p1)]
            founder1.df[idx.p1,] <- alt.alleles[idx.p1]
            founder2.df[idx.p2,] <- alt.alleles[idx.p2]
          } else {
            idx.p1 <- sample(idx.p1.tot, as.numeric(size))
            idx.p2 <- sample(idx.p2.tot, as.numeric(size))
            for(j in 1:length(idx.p1)){
              founder1.df[idx[idx.p1][j],sample(1:ploidy,ploidys[i])] <- alt.alleles[idx[idx.p1][j]]
              founder2.df[idx[idx.p2][j],sample(1:ploidy,ploidys[i])] <- alt.alleles[idx[idx.p2][j]]
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

  command {
    python <<CODE
    from pysam import VariantFile

    bcf_in = VariantFile("~{simulated_vcf}")

    for i in bcf_in.header.samples:
        print(i)
    CODE
  }

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
      --threads=2    
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
    File reads2
    File geno_amb
    File geno_ann
    File geno_bwt
    File geno_pac
    File geno_sa
  }

  command <<<
    /usr/gitc/bwa mem ~{ref} ~{reads1} ~{reads2} | \
      java -jar /usr/gitc/picard.jar SortSam \
        I=/dev/stdin \
        O=~{sampleName}.sorted.bam \
        SORT_ORDER=coordinate \
        CREATE_INDEX=true
    mv ~{sampleName}.sorted.bai ~{sampleName}.sorted.bam.bai
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
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
    File map_file
    Array[File] maternal_trim
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
        
        result <- data.frame(nmk.filt, nmk.id, ok, falso.positivo, ref.ok, alt.ok)
        
        write.table(result, file= paste0(i,".txt"), quote=F, row.names=F, sep="\t")
        
        # tables for mesure depth distribuition
        if(dim(get(i)@gt)[1] != 0){
            idx <- which(strsplit(get(i)@gt[1,1], split=":")[[1]] == "AD")
            if(length(idx) != 0){
            ref.depth <- matrix(sapply(strsplit(sapply(strsplit(get(i)@gt[,-1], split=":"), "[",idx), split=","), "[",1), 
                                ncol = dim(get(i)@gt)[2]-1)
            alt.depth <- matrix(sapply(strsplit(sapply(strsplit(get(i)@gt[,-1], split=":"), "[",idx), split=","), "[",2), 
                                ncol =  dim(get(i)@gt)[2]-1)
            colnames(ref.depth) <- colnames(alt.depth) <- colnames(get(i)@gt[,-1])
            write.table(ref.depth, file = paste0(i,"_ref_depth.txt"), quote=F, row.names=F, sep="\t")
            write.table(alt.depth, file = paste0(i, "_alt_depth.txt"), quote=F, row.names=F, sep="\t")
            } else {
                null.table <- matrix(rep(0,dim(get(i)@gt)[2]-1), ncol=dim(get(i)@gt)[2]-1)
                write.table(null.table, file = paste0(i,"_ref_depth.txt"), quote=F, row.names=F, sep="\t")
                write.table(null.table, file = paste0(i, "_alt_depth.txt"), quote=F, row.names=F, sep="\t")
            }
            # table for GQ 
            idx <- which(strsplit(get(i)@gt[1,1], split=":")$FORMAT == "GQ")
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
        RSCRIPT
    
  >>>

  runtime {
    docker: "taniguti/onemap"
  }

  output {
    File freebayes_pos = "freebayes_site_list.txt"
    File gatk_pos = "gatk_site_list.txt"
    File freebayes_aval_vcf = "freebayes.txt"
    File gatk_aval_vcf = "gatk.txt"
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
      system("cp ~{sep=" "  freebayes_pos} .")
      system("cp ~{sep=" "  gatk_pos} .")
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

task all_maps {

 input {
    File tot_mks
    String methodName
    File simu_vcf
    File vcf_file
    File freebayes_ref_depth2
    File freebayes_alt_depth2
    File gatk_ref_depth2 
    File gatk_alt_depth2
    File gatk_example_alleles
    File freebayes_example_alleles 
    String cross
  }

  command <<<

        R --vanilla --no-save <<RSCRIPT
          library(supermassa4onemap)
          library(onemap)
          library(updog)
          library(reshape2)
          library(vcfR)
          library(doParallel)

          args = commandArgs(trailingOnly=TRUE)

          system("cp ~{sep=" "  freebayes_ref_depth2} .")
          system("cp ~{sep=" "  freebayes_alt_depth2} .")
          system("cp ~{sep=" "  gatk_ref_depth2} .")
          system("cp ~{sep=" "  gatk_alt_depth2} .")
          system("cp ~{gatk_example_alleles} .")
          system("cp ~{freebayes_example_alleles} .")
          method_name <- "~{methodName}"
          tot_mks_file <- "~{tot_mks}"
          simu_vcf_file <- "~{simu_vcf}"
          vcf_file <- "~{vcf_file}"
          cross <- "~{cross}"

          # Functions
          create_filters_report <- function(onemap_obj) {
            segr <- onemap::test_segregation(onemap_obj)
            distorted <- onemap::select_segreg(segr, distorted = T)
            bins <- onemap::find_bins(onemap_obj)
            total_variants <- onemap_obj[[3]]
            filters_tab <- data.frame("n_markers"= total_variants, 
                                      "distorted_markers"=length(distorted), 
                                      "redundant_markers"=total_variants - length(bins))
            return(filters_tab)
          }

          create_maps_report <- function(onemap_obj, tot_mks) {
            assign("onemap_obj", onemap_obj, envir=.GlobalEnv)
            twopts <- rf_2pts(onemap_obj)
            assign("twopts", twopts, envir=.GlobalEnv)
            
            true_mks <- which(onemap_obj[[9]] %in% tot_mks[,2])
            seq_true <- make_seq(twopts, true_mks)
            map_df <- map(seq_true, mds.seq = T)
            while(class(map_df) == "integer"){
              seq_true <- make_seq(twopts, map_df)
              map_df <- map(input.seq = seq_true, mds.seq = T)
            }
            map_info <- data.frame("mk.name"= colnames(onemap_obj[[1]])[map_df[[1]]], 
                                  "pos" = onemap_obj[[9]][map_df[[1]]],
                                  "rf" = c(0,cumsum(haldane(map_df[[3]]))))
            return (map_info)
          }


          create_errors_report <- function(onemap_obj, gab) {
            pos <- which(gab[[9]] %in% onemap_obj[[9]])
            pos.inv <- which(onemap_obj[[9]] %in% gab[[9]])
            gab.pos <- gab[[9]][pos]
            gab.geno <- gab[[1]][,pos]
            colnames(gab.geno) <- gab.pos
            gab.geno <-reshape2::melt(gab.geno)
            colnames(gab.geno) <- c("MK", "POS", "gabGT")
            meth.geno <- onemap_obj[[1]][,pos.inv]
            meth.error <- onemap_obj[[11]][pos.inv + rep(c(0:(onemap_obj[[2]]-1))*onemap_obj[[3]], each=length(pos.inv)),]
            meth.pos <- onemap_obj[[9]][pos.inv]
            colnames(meth.geno) <- meth.pos
            meth.geno <- reshape2::melt(meth.geno)
            colnames(meth.geno) <- c("MK", "POS", "methGT")
            pos.error <- sapply(strsplit(rownames(meth.error), split = "_"), "[",2)
            ind.error <- paste0(sapply(strsplit(rownames(meth.error), split = "_"), "[", 3), "_", sapply(strsplit(rownames(meth.error), split = "_"), "[", 4))
            meth.error <- as.data.frame(cbind(ind.error, pos.error, meth.error))
            colnames(meth.error) <- c("MK", "POS", "A", "AB", "BA", "B")
            error.info <- merge(gab.geno, meth.geno)
            error.info <- merge(error.info, meth.error)
            return (error.info)
          }

          write_report <- function(filters_tab, out_name) {
            write.table(filters_tab, file=out_name, row.names=F, quote=F)
          }


          make_vcf <- function(vcf.old, depths, method){
            # The input od polyRAD need to be a VCF, then this part takes the allele depth from "depths" and put at AD field of input vcf
            idx <- system(paste0("grep -in 'CHROM' ", vcf.old), intern = T) # This part only works in linux OS
            idx.i <- strsplit(idx, split = ":")[[1]][1]
            seed <- sample(1:10000, 1)
            system(paste0("head -n ", idx.i," ", vcf.old, " > head.",seed))
            
            vcf.tab <- read.table(vcf.old, stringsAsFactors = F)
            
            if(all(rownames(depths[[1]]) == paste0(vcf.tab[,1], "_", vcf.tab[,2]))){
              
              vcf.init <- vcf.tab[,1:8]
              AD.colum <- rep("AD", dim(vcf.init)[1])
              vcf.init <- cbind(vcf.init, AD.colum)
              
              rs <- rownames(depths[[1]])
              vcf.init[,3] <- rs
            } else {
              temp.tab <- read.table(paste0(method,"_example4ref_alt_alleles.txt"))
              vcf.init <- cbind(temp.tab[,1:2],paste0(temp.tab[,1], "_", temp.tab[,2]), temp.tab[,3:4], 
                                rep(".", dim(temp.tab)[1]),rep(".", dim(temp.tab)[1]), rep(".", dim(temp.tab)[1]), rep("AD",dim(temp.tab)[1]))
            }
            
            ind.n <- colnames(depths[[1]]) # The names came in different order
            
            header <- strsplit(idx, split = "\t")[[1]]
            ind.vcf <- header[10:length(header)]
            ind.n <- factor(ind.n, levels = ind.vcf)
            
            depths[[1]] <- depths[[1]][,order(ind.n)]
            depths[[2]] <- depths[[2]][,order(ind.n)]
            
            comb.depth <- matrix(paste0(as.matrix(depths[[1]]), ",", as.matrix(depths[[2]])), ncol = ncol(depths[[2]]))
            colnames(comb.depth) <- ind.vcf
            #hmc.file <- cbind(rs, comb.depth)
            
            vcf.body <- cbind(vcf.init, comb.depth)
            
            write.table(vcf.body, file = paste0("temp.body.", seed), quote = FALSE, sep = "\t", row.names = FALSE, col.names = F) 
            
            system(paste0("cat head.",seed," temp.body.",seed," > temp.",seed,".vcf"))
            return(paste0("temp.",seed, ".vcf"))
          }


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
          out_name <- paste0(method_name, "_filters_dfAndGQ.txt")
          filters_tab <- create_filters_report(df)
          write_report(filters_tab, out_name)

          ## MAPS REPORT - DF
          out_name <- paste0(method_name, "_map_df.txt")
          maps_tab <- create_maps_report(df, tot_mks)
          write_report(maps_tab, out_name)

          # MAPS REPORT - GQ
          aval.gq <- extract_depth(vcfR.object=vcf,
                                  onemap.object=df,
                                  vcf.par="GQ",
                                  parent1="P1",
                                  parent2="P2",
                                  f1 = f1,
                                  recovering=FALSE)

          aval.gq <- create_probs(df, genotypes_errors=aval.gq)

          out_name <- paste0(method_name, "_map_GQ.txt")
          maps_gq_tab <- create_maps_report(aval.gq, tot_mks)
          write_report(maps_gq_tab, out_name)

          out_name <- paste0(method_name, "_error_GQ.txt")
          errors_tab <- create_errors_report(aval.gq, gab)
          write_report(errors_tab, out_name)


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
            out_name <- paste0(method_name, "_filters_", metodology, ".txt")
            filters_tab <- create_filters_report(error_aval)
            write_report(filters_tab, out_name)
            
            ## Maps 
            out_name <- paste0(method_name, "_map_", metodology, ".txt")
            maps_tab <- create_maps_report(onemap_obj = error_aval, tot_mks = tot_mks)
            write_report(maps_tab, out_name)
            
            ## Errors info
            out_name <- paste0(method_name, "_error_", metodology, ".txt")
            errors_tab <- create_errors_report(error_aval, gab)
            write_report(errors_tab, out_name)
          }

    
          ## Depths from bam 
          depths.alt <- read.table(paste0(method_name, "_alt_depth_bam.txt"), header = T)
          depths.ref <- read.table(paste0(method_name, "_ref_depth_bam.txt"), header = T)

          depths <- list("ref"=depths.ref, "alt"=depths.alt)

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
            out_name <- paste0(method_name, "_filters_bam_", metodology, ".txt")
            filters_tab <- create_filters_report(error_aval)
            write_report(filters_tab, out_name)
            
            ## Maps 
            out_name <- paste0(method_name, "_map_bam_", metodology, ".txt")
            maps_tab <- create_maps_report(error_aval, tot_mks)
            write_report(maps_tab, out_name)
            
            ## Errors info
            out_name <- paste0(method_name, "_error_bam_", metodology, ".txt")
            errors_tab <- create_errors_report(onemap_obj = error_aval, gab = gab)
            write_report(errors_tab, out_name)
          }
        
        RSCRIPT
  >>>

  runtime{
    docker:"taniguti/onemap"
  }

  output{
    File filters_dfAndGQ = "~{methodName}_filters_dfAndGQ.txt"
    File filters_polyrad = "~{methodName}_filters_polyrad.txt"
    File filters_supermassa = "~{methodName}_filters_supermassa.txt"
    File filters_updog = "~{methodName}_filters_updog.txt"
    File filters_bam_polyrad = "~{methodName}_filters_bam_polyrad.txt"
    File filters_bam_supermassa = "~{methodName}_filters_bam_supermassa.txt"
    File filters_bam_updog = "~{methodName}_filters_bam_updog.txt"
    File map_df = "~{methodName}_map_df.txt"
    File map_GQ = "~{methodName}_map_GQ.txt"
    File map_polyrad = "~{methodName}_map_polyrad.txt"
    File map_supermassa = "~{methodName}_map_supermassa.txt"
    File map_updog = "~{methodName}_map_updog.txt"
    File map_bam_polyrad = "~{methodName}_map_bam_polyrad.txt"
    File map_bam_supermassa = "~{methodName}_map_bam_supermassa.txt"
    File map_bam_updog = "~{methodName}_map_bam_updog.txt"
    File error_info_GQ = "~{methodName}_error_GQ.txt"
    File error_info_updog = "~{methodName}_error_updog.txt"
    File error_info_polyrad = "~{methodName}_error_polyrad.txt"
    File error_info_supermassa = "~{methodName}_error_supermassa.txt"
    File error_info_bam_updog = "~{methodName}_error_bam_updog.txt"
    File error_info_bam_polyrad = "~{methodName}_error_bam_polyrad.txt"
    File error_info_bam_supermassa = "~{methodName}_error_bam_supermassa.txt"

  }
}

