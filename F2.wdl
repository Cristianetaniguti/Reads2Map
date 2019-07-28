version 1.0

import "./structs/f2.wdl"

workflow SimulateF2 {

  input {
    ReferenceFasta references
    Family family
  }
    
  call GenerateAlternativeGenome {
    input:
    seed=family.seed,
    ref_genome = references.ref_fasta
  }

  call CreatePedigreeSimulatorInputs {
    input:
      seed = family.seed,
      snps = GenerateAlternativeGenome.snps,
      indels = GenerateAlternativeGenome.indels,
      cmBymb = family.cmBymb,
      ref = references.ref_fasta,
      ref_fai = references.ref_fasta_index
  }

  call RunPedigreeSimulator {
    input:
      map_file = CreatePedigreeSimulatorInputs.mapfile,
      founder_file = CreatePedigreeSimulatorInputs.founderfile,
      chrom_file = CreatePedigreeSimulatorInputs.chromfile,
      par_file = CreatePedigreeSimulatorInputs.parfile
  }

  call PedgreeSimulator2vcf {
    input:
      genotypes_dat = RunPedigreeSimulator.genotypes_dat,
      map_file = CreatePedigreeSimulatorInputs.mapfile,
      chrom_file = CreatePedigreeSimulatorInputs.chromfile,
      tot_mks = CreatePedigreeSimulatorInputs.tot_mks
  }

  call GenerateSampleNames {
    input: 
      samples = family.samples
  }

  scatter (sampleName in GenerateSampleNames.names) {

    call RunVcf2diploid {
      input:
        sampleName = sampleName,
        ref_genome = references.ref_fasta,
        simu_vcf = PedgreeSimulator2vcf.simu_vcf
    }

    call SimulateRADseq {
      input:
        enzyme = family.enzyme,
        sampleName = sampleName,
        maternal_genomes = RunVcf2diploid.maternal_genomes,
        paternal_genomes = RunVcf2diploid.paternal_genomes
    }

    call SimulateIlluminaReads {
      input:
        maternal_trim = SimulateRADseq.maternal_trim,
        paternal_trim = SimulateRADseq.paternal_trim,
        sampleName = sampleName
    }

    call RunBwaAlignment {
      input:
        sampleName = sampleName,
        reads1 = SimulateIlluminaReads.reads1,
        reads2 = SimulateIlluminaReads.reads2,
        ref = references.ref_fasta,
        geno_amb = references.ref_amb,
        geno_ann = references.ref_ann,
        geno_bwt = references.ref_bwt,
        geno_pac = references.ref_pac,
        geno_sa = references.ref_sa
    }

    call AddAlignmentHeader {
      input:
        sampleName = sampleName,
        bam_file = RunBwaAlignment.bam_file,
        bam_idx = RunBwaAlignment.bam_idx
    }

    call HaplotypeCallerERC {
      input:
        ref = references.ref_fasta,
        geno_fai = references.ref_fasta_index,
        sampleName = sampleName,
        bam_rg = AddAlignmentHeader.bam_rg,
        bam_rg_idx = AddAlignmentHeader.bam_rg_index,
        geno_dict = references.ref_dict
    }
  }

  call CreateGatkDatabase {
    input:
      path_gatkDatabase = "my_database",
      GVCFs = HaplotypeCallerERC.GVCF,
      GVCFs_idx = HaplotypeCallerERC.GVCF_idx
  }

  call GenotypeGVCFs {
    input:
      workspace_tar = CreateGatkDatabase.workspace_tar,
      output_vcf_filename = family.name + "_gatk.vcf",
      ref = references.ref_fasta,
      geno_fai = references.ref_fasta_index,
      geno_dict = references.ref_dict
  }

  call RunFreebayes {
    input:
      freebayesVCFname = family.name + "_freebayes.vcf",
      ref = references.ref_fasta,
      bam_rg = AddAlignmentHeader.bam_rg
  }

  call VcftoolsApplyFilters{
    input:
      freebayesVCF = RunFreebayes.freebayesVCF,
      gatkVCF = GenotypeGVCFs.gatkVCF
  }

  call CalculateVcfMetrics {
    input:
      freebayesVCF = VcftoolsApplyFilters.freebayesVCF_F,
      gatkVCF = VcftoolsApplyFilters.gatkVCF_F,
      tot_mks = CreatePedigreeSimulatorInputs.tot_mks,
      map_file = CreatePedigreeSimulatorInputs.mapfile,
      maternal_trim = SimulateRADseq.maternal_trim
  }

  Array[String] methods = ["gatk", "freebayes"]
  Array[File] vcfs = [VcftoolsApplyFilters.gatkVCF_F, VcftoolsApplyFilters.freebayesVCF_F]
  Array[Pair[String, File]] program_and_vcf = zip(methods, vcfs)

  scatter (vcf in program_and_vcf){
    call all_maps{
      input:
        tot_mks = CreatePedigreeSimulatorInputs.tot_mks,
        simu_vcf = PedgreeSimulator2vcf.simu_vcf,
        methodName = vcf.left,
        vcf_file = vcf.right
    }
  }

  output {
    File freebayes_vcf = VcftoolsApplyFilters.freebayesVCF_F
    File gatk_vcf = VcftoolsApplyFilters.gatkVCF_F
    File freebayes_aval_vcf = CalculateVcfMetrics.freebayes_aval_vcf
    File gatk_aval_vcf = CalculateVcfMetrics.gatk_aval_vcf
    File freebayes_ref_depth = CalculateVcfMetrics.freebayes_ref_depth
    File freebayes_alt_depth = CalculateVcfMetrics.freebayes_alt_depth
    File gatk_ref_depth = CalculateVcfMetrics.gatk_ref_depth
    File gatk_alt_depth = CalculateVcfMetrics.gatk_alt_depth
    File tot_mks = CreatePedigreeSimulatorInputs.tot_mks
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
  }

  command <<<
    python /opt/scripts/pedsim_files.py --indels ~{indels} --snps ~{snps} --reference ~{ref} --seed ~{seed}
  >>>

  runtime {
    docker: "taniguti/miniconda-alpine"
  }

  output {
    File mapfile = "markers.txt"
    File founderfile = "founders.txt"
    File parfile = "parameters.txt"
    File chromfile = "chromosome.txt"
    File tot_mks = "tot_mks.txt"
  }

}

task RunPedigreeSimulator {
  input {
    File map_file
    File founder_file
    File chrom_file
    File par_file
  }

  command <<<
    sed -i 's+inb.chrom+~{chrom_file}+g' ~{par_file}
    sed -i 's+mapfile.map+~{map_file}+g' ~{par_file}
    sed -i 's+founderfile.gen+~{founder_file}+g' ~{par_file}
    java -jar /usr/jars/PedigreeSim.jar ~{par_file}
  >>>

  runtime {
    docker: "taniguti/java-in-the-cloud"
  }

  output {
    File genotypes_dat = "sim_inb_genotypes.dat"
  }

}

# Parse pedsim output (.dat) into VCF
task PedgreeSimulator2vcf {

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
    print(mks)
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
    Int samples
  }

  command {
    python <<CODE
    samples = ~{samples}

    names = ["P1", "P2", "F1"] + ["F2_%03d" % i for i in range(1, samples + 1)]
    for i in names:
        print(i)
    CODE
  }

  runtime {
    docker: "python:3.7"
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
    String sampleName
  }

  command <<<
    set -e
    /pirs/src/pirs/pirs simulate \
      --diploid ~{maternal_trim} ~{paternal_trim} \
      --read-len=100 \
      --coverage=100 \
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
      -O ~{sampleName}_rawLikelihoods.g.vcf
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
    File freebayes_aval_vcf = "freebayes.txt"
    File gatk_aval_vcf = "gatk.txt"
    File freebayes_ref_depth = "freebayes_ref_depth.txt"
    File freebayes_alt_depth = "freebayes_alt_depth.txt"
    File gatk_ref_depth = "gatk_ref_depth.txt"
    File gatk_alt_depth = "gatk_alt_depth.txt"
  }
}

task all_maps {

 input {
    File tot_mks
    String methodName
    File simu_vcf
    File vcf_file
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

          method_name <- "~{methodName}"
          tot_mks_file <- "~{tot_mks}"
          simu_vcf_file <- "~{simu_vcf}"
          vcf_file <- "~{vcf_file}"

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
              map_df <- map(seq_true)
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
              meth.error <- cbind(rownames(meth.error), meth.error)
              error.info <- merge(gab.geno, meth.geno)
              error.info <- error.info[order(error.info[,1], error.info[,2]),]
              error.info <- cbind(error.info, meth.error)
              return (error.info)
          }

          write_report <- function(filters_tab, out_name) {
              write.table(filters_tab, file=out_name, row.names=F, quote=F)
          }

          ## KNOWN VARIANTS
          tot_mks <- read.table(tot_mks_file)

          # READING DATA FROM SIMULATED POPULATION
          simu <- read.vcfR(simu_vcf_file)
          gab <- onemap_read_vcfR(vcfR.object=simu,
                                  cross="f2 intercross",
                                  parent1="P1",
                                  parent2="P2",
                                  f1="F1")

          ## READING FINAL VCF FROM PIPELINE
          vcf <- read.vcfR(vcf_file)
          df <- onemap_read_vcfR(vcfR.object=vcf, 
                                cross="f2 intercross", 
                                parent1="P1", 
                                parent2="P2", 
                                f1="F1")

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
                                  f1="F1",
                                  recovering=FALSE)

          aval.gq <- create_probs(df, genotypes_errors=aval.gq)

          out_name <- paste0(method_name, "_map_GQ.txt")
          maps_gq_tab <- create_maps_report(aval.gq, tot_mks)
          write_report(maps_gq_tab, out_name)

          out_name <- paste0(method_name, "_error_info_GQ.txt")
          errors_tab <- create_errors_report(aval.gq, gab)
          write_report(errors_tab, out_name)


          # OTHER TOOLS
          updog.aval <- updog_error(
              vcfR.object=vcf,
              onemap.object=df,
              vcf.par="AD",
              parent1="P1",
              parent2="P2",
              f1="F1",
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
              f1="F1",
              recovering = TRUE,
              mean_phred = 20,
              cores = 3,
              depths = NULL)

          polyrad.aval <- polyRAD_error(
              vcf=vcf_file, 
              onemap.obj=df,
              parent1="P1",
              parent2="P2",
              f1="F1",
              crosstype="f2 intercross")

          metodologies <- list(updog = updog.aval, supermassa = supermassa.aval, polyrad = polyrad.aval)
          for (metodology in names(metodologies)){
              error_aval <- metodologies[[metodology]]
              ## Filters
              out_name <- paste0(method_name, "_filters_", metodology, ".txt")
              filters_tab <- create_filters_report(error_aval)
              write_report(filters_tab, out_name)

              ## Maps updog
              out_name <- paste0(method_name, "_map_", metodology, ".txt")
              maps_tab <- create_maps_report(error_aval, tot_mks)
              write_report(maps_tab, out_name)

              ## Errors info
              out_name <- paste0(method_name, "_error_", metodology, ".txt")
              errors_tab <- create_errors_report(error_aval, gab)
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
    File map_df = "~{methodName}_map_df.txt"
    File map_GQ = "~{methodName}_map_GQ.txt"
    File map_polyrad = "~{methodName}_map_polyrad.txt"
    File map_supermassa = "~{methodName}_map_supermassa.txt"
    File error_info_GQ = "~{methodName}_error_info_GQ.txt"
    File error_info_updog = "~{methodName}_error_updog.txt"
    File error_info_polyrad = "~{methodName}_error_polyrad.txt"
    File error_info_supermassa = "~{methodName}_error_supermassa.txt"
  }
}
