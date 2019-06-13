version 1.0

import "./structs/f2.wdl"

workflow F2 {
  input {
    ReferenceFasta references
    Family family
  }

  String path_gatkDatabase = "my_database"
  Array[String] sampleNames = read_lines(family.samples_names_file)
    
  call create_alt_genome {
    input:
    seed=family.seed,
    ref_genome = references.ref_fasta
  }

  call pedsim_files {
    input:
      seed=family.seed,
      snps = create_alt_genome.snps,
      indels = create_alt_genome.indels,
      cmBymb = family.cmBymb,
      ref = references.ref_fasta,
      ref_fai = references.ref_fasta_index
  }

  call pedigreeSim {
    input:
      map_file = pedsim_files.mapfile,
      founder_file = pedsim_files.founderfile,
      chrom_file = pedsim_files.chromfile,
      par_file = pedsim_files.parfile
  }

  call pedsim2vcf {
    input:
      genotypes_dat = pedigreeSim.genotypes_dat,
      map_file = pedsim_files.mapfile,
      chrom_file = pedsim_files.chromfile,
      tot_mks = pedsim_files.tot_mks
  }


  scatter (sampleName in sampleNames) {
    call reads_simulations {
      input:
        seed=family.seed,
        maternal_trim = create_frags.maternal_trim,
        paternal_trim = create_frags.paternal_trim,
        sampleName = sampleName
    }

    call vcf2diploid {
      input:
        sampleName = sampleName,
        ref_genome = references.ref_fasta,
        simu_vcf = pedsim2vcf.simu_vcf
    }

    call create_frags {
      input:
        enzyme = family.enzyme,
        sampleName = sampleName,
        maternal_genomes = vcf2diploid.maternal_genomes,
        paternal_genomes = vcf2diploid.paternal_genomes
    }

    call alignment {
      input:
        sampleName = sampleName,
        reads1 = reads_simulations.reads1,
        reads2 = reads_simulations.reads2,
        ref = references.ref_fasta,
        geno_amb = references.ref_amb,
        geno_ann = references.ref_ann,
        geno_bwt = references.ref_bwt,
        geno_pac = references.ref_pac,
        geno_sa = references.ref_sa
    }

    call add_labs {
      input:
        sampleName = sampleName,
        bam_file = alignment.bam_file,
        bam_idx = alignment.bam_idx
    }

    call HaplotypeCallerERC {
      input:
        ref = references.ref_fasta,
        geno_fai = references.ref_fasta_index,
        sampleName = sampleName,
        bam_rg = add_labs.bam_rg,
        bam_rg_idx = add_labs.bam_rg_index,
        geno_dict = references.ref_dict
    }
  }

  call GenotypeGVCFs {
    input:
      workspace_tar = create_gatk_database.workspace_tar,
      output_vcf_filename = family.name + "_gatk.vcf",
      ref = references.ref_fasta,
      geno_fai = references.ref_fasta_index,
      geno_dict = references.ref_dict
  }

  call create_gatk_database {
    input:
      path_gatkDatabase = path_gatkDatabase,
      GVCFs = HaplotypeCallerERC.GVCF,
      GVCFs_idx = HaplotypeCallerERC.GVCF_idx
  }

  call freebayes {
    input:
      freebayesVCFname = family.name + "_freebayes.vcf",
      ref = references.ref_fasta,
      bam_rg = add_labs.bam_rg
  }

  call create_popmapFile {
    input:
      sampleNamesFile = family.samples_names_file
  }

  call ref_map {
    input:
      bam_rg = add_labs.bam_rg,
      popmapfile = create_popmapFile.popmapfile
  }

  call aval_vcf {
    input:
      stacksVCF = ref_map.stacksVCF,
      freebayesVCF = freebayes.freebayesVCF,
      gatkVCF = GenotypeGVCFs.gatkVCF,
      tot_mks = pedsim_files.tot_mks,
      map_file = pedsim_files.mapfile,
      maternal_trim = create_frags.maternal_trim
  }

  output {
    File stacks_vcf = ref_map.stacksVCF
    File freebayes_vcf = freebayes.freebayesVCF
    File gatk_vcf = GenotypeGVCFs.gatkVCF
    File freebayes_aval_vcf = aval_vcf.freebayes_aval_vcf
    File gatk_aval_vcf = aval_vcf.gatk_aval_vcf
    File stacks_aval_vcf = aval_vcf.stacks_aval_vcf
    File freebayes_ref_depth = aval_vcf.freebayes_ref_depth
    File freebayes_alt_depth = aval_vcf.freebayes_alt_depth
    File gatk_ref_depth = aval_vcf.gatk_ref_depth
    File gatk_alt_depth = aval_vcf.gatk_alt_depth
    File stacks_ref_depth = aval_vcf.stacks_ref_depth
    File stacks_alt_depth = aval_vcf.stacks_alt_depth
    File vcfRs = aval_vcf.vcfRs
    File tot_mks = pedsim_files.tot_mks
  }
}

# Creates homologous genome with some variation
# specified with -s and -d
task create_alt_genome {
  input {
    Int seed
    File ref_genome
  }

  output {
    File alt_fasta = "alt.snp.indel.fa"
    File indels = "alt.indel.lst"
    File snps = "alt.snp.lst"
  }
  command <<<

        /pirs/src/pirs/pirs diploid ~{ref_genome} -s 0.001 -d 0.0001 -v 0 -o alt --random-seed ~{seed}
    
  >>>
  runtime {
    docker: "pirs-ddrad-cutadapt:v1"
  }

}

task pedsim_files {
  input {
    File snps
    File indels
    Float cmBymb
    File ref
    File ref_fai
    Int seed
  }

  output {
    File mapfile = "mapfile.map"
    File founderfile = "founderfile.gen"
    File parfile = "sim.par"
    File chromfile = "inb.chrom"
    File tot_mks = "tot_mks.txt"
  }
  command <<<

        R --vanilla --no-save <<RSCRIPT
        snps <- read.table("~{snps}", stringsAsFactors = FALSE)
        indels <- read.table("~{indels}", stringsAsFactors = FALSE)
        pos.ref <- indels[,2]
        sinal <- indels[,4]

        # Nos arquivos de saída do pirs nao consta a ultima base antes do polimorfismo
        # Quanto se trata de indels negativos para a referencia, a posição anterior a apontada
        # é a ultima base antes do polimorfismo
        pos.ref[which(sinal=="-")] <- pos.ref[which(sinal=="-")] -1

        # search last base before the indels (information needed by VCF)
        command  <- c(paste("samtools faidx ~{ref}"),paste0("Chr10:",pos.ref,"-",pos.ref))
        bases <- system(paste0(command, collapse = " "), intern = T)

        bases.bf <- matrix(bases, ncol=2, byrow = T)[,2]
        alt <- bases.bf
        tmp <- paste0(bases.bf[which(sinal=="+")], indels[,6][which(sinal=="+")])
        alt[which(sinal=="+")] <- tmp
        ref <- bases.bf
        tmp <- paste0(bases.bf[which(sinal=="-")], indels[,6][which(sinal=="-")])
        ref[which(sinal=="-")] <- tmp

        # a posição no vcf e no mapa são em relação ao genoma de referência
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
        write.table(map_file, file = paste0("mapfile.map"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

        founder_file <- data.frame(marker=marker, P1_1=tot.mks[,3] , P1_2=tot.mks[,3], P2_1=tot.mks[,4], P2_2=tot.mks[,4])
        write.table(founder_file, file = paste0("founderfile.gen"), quote=FALSE, col.names = TRUE, row.names = FALSE, sep = "\t" )

        ## Parameters file
        parameter <- paste0("PLOIDY = 2
                                        MAPFUNCTION = HALDANE
                                        MISSING = NA
                                        CHROMFILE = inb.chrom
                                        POPTYPE = F2
                                        SEED = ~{seed}
                                        POPSIZE = 150
                                        MAPFILE = mapfile.map
                                        FOUNDERFILE = founderfile.gen
                                        OUTPUT = sim_inb")

        write.table(parameter, file = paste0("sim.par"), quote=FALSE, col.names = FALSE, row.names = FALSE, sep = "\t" )
        chrom <- data.frame("chromosome"= "Chr10", "length"= pos.map[which.max(pos.map)], "centromere"=pos.map[which.max(pos.map)]/2, "prefPairing"= 0.0, "quadrivalents"=0.0)
        write.table(chrom, file= "inb.chrom", quote = F, col.names = T, row.names = F, sep= "\t")
        RSCRIPT
    
  >>>
  runtime {
    docker: "r-samtools:v1"
  }
}

task pedigreeSim {
  input {
    File map_file
    File founder_file
    File chrom_file
    File par_file
  }

  output {
    File genotypes_dat = "sim_inb_genotypes.dat"
  }
  command <<<

        sed -i 's+inb.chrom+~{chrom_file}+g' ~{par_file}
        sed -i 's+mapfile.map+~{map_file}+g' ~{par_file}
        sed -i 's+founderfile.gen+~{founder_file}+g' ~{par_file}

        java -jar /usr/jars/PedigreeSim.jar ~{par_file}
    
  >>>
  runtime {
    docker: "java-in-the-cloud:v1"
  }

}

# Parse pedsim output (.dat) into VCF
task pedsim2vcf {
  input {
    File genotypes_dat
    File map_file
    File chrom_file
    File tot_mks
  }

  output {
    File simu_vcf = "simu.vcf"
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
    docker: "onemap:v1"
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


  output {
    File GVCF = "${sampleName}_rawLikelihoods.g.vcf"
    File GVCF_idx = "${sampleName}_rawLikelihoods.g.vcf.idx"
  }
  command <<<

        /gatk/gatk HaplotypeCaller \
            -ERC GVCF \
            -R ~{ref} \
            -I ~{bam_rg} \
            -O ~{sampleName}_rawLikelihoods.g.vcf
    
  >>>
  runtime {
    docker: "gatk-picard:v1"
  }

}


# Simulates an Illumina sequencing experiment
# TODO: produce single sample reads
# Fix is required because in some cases pirs fails to
# write nucleotide qualities in fastq file
task reads_simulations {
  input {
    File maternal_trim
    File paternal_trim
    String sampleName
    Int seed
  }


  output {
    File reads1 = "${sampleName}_100_150_1_fix.fq"
    File reads2 = "${sampleName}_100_150_2_fix.fq"
  }
  command <<<

        /pirs/src/pirs/pirs simulate \
            --diploid ~{maternal_trim} ~{paternal_trim} \
            -l 100 -x 50 -m 150 -o ~{sampleName} --random-seed ~{seed}

        /cleanFastq/fixFastq "~{sampleName}_100_150_1.fq" "~{sampleName}_100_150_1_fix.fq" 
        /cleanFastq/fixFastq "~{sampleName}_100_150_2.fq" "~{sampleName}_100_150_2_fix.fq"
    
  >>>
  runtime {
    docker: "pirs-ddrad-cutadapt:v1"
  }

}

# Variant calling using freebayes
task freebayes {
  input {
    String freebayesVCFname
    File ref
    Array[File] bam_rg
  }


  output {
    File freebayesVCF = "${freebayesVCFname}"
  }
  command <<<

        freebayes --genotype-qualities -f ~{ref} ~{sep=" "  bam_rg} > ~{freebayesVCFname}
    
  >>>
  runtime {
    docker: "freebayes:v1"
  }

}



# Insert into a fasta sequence the variants present in a VCF file
task vcf2diploid {
  input {
    String sampleName
    File ref_genome
    File simu_vcf
  }


  output {
    File maternal_genomes = "Chr10_${sampleName}_maternal.fa"
    File paternal_genomes = "Chr10_${sampleName}_paternal.fa"
  }
  command <<<

        java -jar /usr/jars/vcf2diploid.jar -id ~{sampleName} -chr ~{ref_genome} -vcf ~{simu_vcf}
    
  >>>
  runtime {
    docker: "java-in-the-cloud:v1"
  }

}

# Add info to alignment header
task add_labs {
  input {
    String sampleName
    File bam_file
    File bam_idx
  }


  output {
    File bam_rg = "${sampleName}_rg.bam"
    File bam_rg_index = "${sampleName}_rg.bam.bai"
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
    docker: "gatk-picard:v1"
  }
}

# Alignment using bwa mem
task alignment {
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


  output {
    File bam_file = "${sampleName}.sorted.bam"
    File bam_idx = "${sampleName}.sorted.bam.bai"
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

  output {
    File gatkVCF = "${output_vcf_filename}"
    File gatkVCF_index = "${output_vcf_filename}.idx"
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
    docker: "gatk-picard:v1"
  }

}

task aval_vcf {
  input {
    File stacksVCF
    File freebayesVCF
    File gatkVCF
    File tot_mks
    File map_file
    Array[File] maternal_trim
  }

  output {
    File freebayes_aval_vcf = "freebayes.txt"
    File gatk_aval_vcf = "gatk.txt"
    File stacks_aval_vcf = "stacks.txt"
    File freebayes_ref_depth = "freebayes_ref_depth.txt"
    File freebayes_alt_depth = "freebayes_alt_depth.txt"
    File gatk_ref_depth = "gatk_ref_depth.txt"
    File gatk_alt_depth = "gatk_alt_depth.txt"
    File stacks_ref_depth = "stacks_ref_depth.txt"
    File stacks_alt_depth = "stacks_alt_depth.txt"
    File vcfRs = "vcfRs.RData"
  }
  command <<<

        R --vanilla --no-save <<RSCRIPT

        system("cp ~{sep=" "  maternal_trim} .")
        library(vcfR)
        freebayes <- read.vcfR("~{freebayesVCF}")
        gatk <- read.vcfR("~{gatkVCF}")
        stacks <- read.vcfR("~{stacksVCF}")
        vcfRs <- list("gatk"=gatk, "freebayes"=freebayes, "stacks"=stacks)
        save(vcfRs, file="vcfRs.RData")
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

        methods <- c("freebayes", "gatk", "stacks")
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
    docker: "onemap:v1"
  }

}

# SNP calling program
task ref_map {
  input {
    Array[File] bam_rg
    File popmapfile
  }


  output {
    File stacksVCF = "populations.snps.vcf"
  }
  command <<<

        cp ~{sep=" "  bam_rg} .
        ref_map.pl --samples . --popmap ~{popmapfile} -o . -X "populations:--ordered-export --vcf "
    
  >>>
  runtime {
    docker: "stacks:v1"
  }

}

task create_gatk_database {
  input {
    String path_gatkDatabase
    Array[File] GVCFs
    Array[File] GVCFs_idx
  }


  output {
    File workspace_tar = "${path_gatkDatabase}.tar"
  }
  command <<<

        /gatk/gatk GenomicsDBImport \
            --genomicsdb-workspace-path ~{path_gatkDatabase} \
            -L Chr10 \
            -V ~{sep=" -V "  GVCFs} 

        tar -cf ~{path_gatkDatabase}.tar ~{path_gatkDatabase}
    
  >>>
  runtime {
    docker: "gatk-picard:v1"
  }

}

# Simulates RADseq experiment where certain enzyme is
# used to fragment the sequence and then the first X
# bases of each resulting fragment is sequenced.
task create_frags {
  input {
    String enzyme
    String sampleName
    File maternal_genomes
    File paternal_genomes
  }

  output {
    File maternal_frags = "${sampleName}_maternal_fragments.fasta"
    File paternal_frags = "${sampleName}_paternal_fragments.fasta"
    File maternal_stats = "${sampleName}_maternal_statistics.txt"
    File paternal_stats = "${sampleName}_paternal_statistics.txt"
    File maternal_trim = "${sampleName}_maternal_trim.fa"
    File paternal_trim = "${sampleName}_paternal_trim.fa"
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
    docker: "pirs-ddrad-cutadapt:v1"
  }

}


# Creating input files
task create_popmapFile {
  input {
    File sampleNamesFile
  }

  output {
    File popmapfile = "popmap.txt"
  }

  command <<<

        R --vanilla --no-save <<RSCRIPT
        names <- read.table("~{sampleNamesFile}", header = F)
        mapnames <- paste0(t(names), '_rg')
        mapdf <- data.frame(mapnames, rep(1, length(mapnames)))
        write.table(mapdf, file = 'popmap.txt', sep = '\t', col.names = F, row.names = F, quote=F)
        RSCRIPT
    
  >>>
  runtime {
    docker: "onemap:v1"
  }

}
