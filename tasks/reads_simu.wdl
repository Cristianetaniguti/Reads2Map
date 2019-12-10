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
        real_phases               = CreatePedigreeSimulatorInputs.real_phases
    }
  }

  call CreateTables{
    input:
        times                     = CreateMaps.times,
        cmBymb                    = family.cmBymb,
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
        map_df                    = CreateMaps.map_df,
        map_GQ                    = CreateMaps.map_GQ,
        map_polyrad               = CreateMaps.map_polyrad,
        map_supermassa            = CreateMaps.map_supermassa,
        map_updog                 = CreateMaps.map_updog,
        map_gusmap                = CreateMaps.map_gusmap,
        map_bam_polyrad           = CreateMaps.map_bam_polyrad,
        map_bam_supermassa        = CreateMaps.map_bam_supermassa,
        map_bam_updog             = CreateMaps.map_bam_updog,
        map_bam_gusmap            = CreateMaps.map_bam_gusmap,
        error_info_df             = CreateMaps.error_info_df,
        error_info_GQ             = CreateMaps.error_info_GQ,
        error_info_updog          = CreateMaps.error_info_updog,
        error_info_polyrad        = CreateMaps.error_info_polyrad,
        error_info_supermassa     = CreateMaps.error_info_supermassa,
        error_info_bam_updog      = CreateMaps.error_info_bam_updog,
        error_info_bam_polyrad    = CreateMaps.error_info_bam_polyrad,
        error_info_bam_supermassa = CreateMaps.error_info_bam_supermassa,
        filters_dfAndGQ           = CreateMaps.filters_dfAndGQ,
        filters_polyrad           = CreateMaps.filters_polyrad,
        filters_supermassa        = CreateMaps.filters_supermassa,
        filters_updog             = CreateMaps.filters_updog,
        filters_bam_polyrad       = CreateMaps.filters_bam_polyrad,
        filters_bam_supermassa    = CreateMaps.filters_bam_supermassa,
        filters_bam_updog         = CreateMaps.filters_bam_updog
    }

  output {
    File data1_depths_geno_prob   = CreateTables.data1_depths_geno_prob
    File data2_maps               = CreateTables.data2_maps
    File data3_coverage           = CreateTables.data3_coverage
    File data4_filters            = CreateTables.data4_filters
    File data5_SNPcall_efficiency = CalculateVcfMetrics.data5_SNPcall_efficiency
    File data6_times              = CreateTables.data6_times
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
          real_phases <- read.table("~{real_phases}")
          source("/opt/scripts/functions.R")

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
          write_report(filters_tab[[1]], out_name)

          ## MAPS REPORT - DF
          out_name <- paste0(method_name, "_map_df.txt")
          times <-system.time(maps_tab <- create_maps_report(filters_tab[[2]], tot_mks, gab))
          write_report(maps_tab, out_name)
          times <- data.frame(meth =paste0(method_name, "_map_df"), time = times[3])

          out_name <- paste0(method_name, "_error_df.txt")
          errors_tab <- create_errors_report(df, gab)
          write_report(errors_tab, out_name)

          # MAPS REPORT - GQ
          aval.gq <- extract_depth(vcfR.object=vcf,
                                  onemap.object=df,
                                  vcf.par="GQ",
                                  parent1="P1",
                                  parent2="P2",
                                  f1 = f1,
                                  recovering=FALSE)

          aval.gq <- create_probs(df, genotypes_errors=aval.gq)

          filters_tab <- create_filters_report(aval.gq)
          out_name <- paste0(method_name, "_map_GQ.txt")
          times_temp <- system.time(maps_gq_tab <- create_maps_report(filters_tab[[2]], tot_mks, gab))
          write_report(maps_gq_tab, out_name)
          times_temp <- data.frame(meth =paste0(method_name, "_map_GQ"), time = times_temp[3])
          times <- rbind(times, times_temp)

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
            write_report(filters_tab[[1]], out_name)

            ## Maps
            out_name <- paste0(method_name, "_map_", metodology, ".txt")
            times_temp <- system.time(maps_tab <- create_maps_report(input.seq = filters_tab[[2]], tot_mks = tot_mks, gab))
            write_report(maps_tab, out_name)
            times_temp <- data.frame(meth =paste0(method_name, "_map_",metodology), time = times_temp[3])
            times <- rbind(times, times_temp)

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
            write_report(filters_tab[[1]], out_name)

            ## Maps
            out_name <- paste0(method_name, "_map_bam_", metodology, ".txt")
            times_temp <- system.time(maps_tab <- create_maps_report(filters_tab[[2]], tot_mks = tot_mks, gab))
            write_report(maps_tab, out_name)
            times_temp <- data.frame(meth =paste0(method_name, "_map_bam_",metodology), time = times_temp[3])
            times <- rbind(times, times_temp)

            ## Errors info
            out_name <- paste0(method_name, "_error_bam_", metodology, ".txt")
            errors_tab <- create_errors_report(onemap_obj = error_aval, gab = gab)
            write_report(errors_tab, out_name)
          }

          ## Gusmap maps
          out_name <- paste0(method_name, "_map_gusmap.txt")
          times_temp <- system.time(map_gus <- create_gusmap_report(vcf_file, gab))
          write_report(map_gus, out_name)
          times_temp <- data.frame(meth =paste0(method_name, "_map_gusmap"), time = times_temp[3])
          times <- rbind(times, times_temp)

          out_name <- paste0(method_name, "_map_bam_gusmap.txt")
          times_temp <- system.time(map_gus <- create_gusmap_report(new.vcf, gab))
          write_report(map_gus, out_name)
          times_temp <- data.frame(meth =paste0(method_name, "_map_bam_gusmap"), time = times_temp[3])
          times <- rbind(times, times_temp)
          write.table(times, paste0(method_name,"_times.txt"))

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
    File map_gusmap = "~{methodName}_map_gusmap.txt"
    File map_bam_gusmap = "~{methodName}_map_bam_gusmap.txt"
    File error_info_df = "~{methodName}_error_df.txt"
    File error_info_GQ = "~{methodName}_error_GQ.txt"
    File error_info_updog = "~{methodName}_error_updog.txt"
    File error_info_polyrad = "~{methodName}_error_polyrad.txt"
    File error_info_supermassa = "~{methodName}_error_supermassa.txt"
    File error_info_bam_updog = "~{methodName}_error_bam_updog.txt"
    File error_info_bam_polyrad = "~{methodName}_error_bam_polyrad.txt"
    File error_info_bam_supermassa = "~{methodName}_error_bam_supermassa.txt"
    File times = "~{methodName}_times.txt"

  }
}

task CreateTables{
  input{
    Array[File] times
    Float cmBymb
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
    Array[File] map_df
    Array[File] map_GQ
    Array[File] map_polyrad
    Array[File] map_supermassa
    Array[File] map_updog
    Array[File] map_bam_polyrad
    Array[File] map_bam_supermassa
    Array[File] map_bam_updog
    Array[File] map_gusmap
    Array[File] map_bam_gusmap
    Array[File] error_info_df
    Array[File] error_info_GQ
    Array[File] error_info_updog
    Array[File] error_info_polyrad
    Array[File] error_info_supermassa
    Array[File] error_info_bam_updog
    Array[File] error_info_bam_polyrad
    Array[File] error_info_bam_supermassa
    Array[File] filters_dfAndGQ
    Array[File] filters_polyrad
    Array[File] filters_supermassa
    Array[File] filters_updog
    Array[File] filters_bam_polyrad
    Array[File] filters_bam_supermassa
    Array[File] filters_bam_updog
  }

  command <<<

        R --vanilla --no-save <<RSCRIPT

          library(reshape2)
          library(vcfR)
          system("cp ~{sep= " " times } .")
          system("cp ~{sep= " " map_df } .")
          system("cp ~{sep= " " map_GQ } .")
          system("cp ~{sep= " " map_polyrad } .")
          system("cp ~{sep= " " map_supermassa } .")
          system("cp ~{sep= " " map_updog } .")
          system("cp ~{sep= " " map_bam_polyrad } .")
          system("cp ~{sep= " " map_bam_supermassa } .")
          system("cp ~{sep= " " map_bam_updog } .")
          system("cp ~{sep= " " map_gusmap } .")
          system("cp ~{sep= " " map_bam_gusmap } .")
          system("cp ~{sep= " " error_info_df } .")
          system("cp ~{sep= " " error_info_GQ } .")
          system("cp ~{sep= " " error_info_updog } .")
          system("cp ~{sep= " " error_info_polyrad } .")
          system("cp ~{sep= " " error_info_supermassa } .")
          system("cp ~{sep= " " error_info_bam_updog } .")
          system("cp ~{sep= " " error_info_bam_polyrad } .")
          system("cp ~{sep= " " error_info_bam_supermassa } .")
          system("cp ~{sep = " " filters_dfAndGQ } .")
          system("cp ~{sep = " " filters_polyrad } .")
          system("cp ~{sep = " " filters_supermassa } .")
          system("cp ~{sep = " " filters_updog } .")
          system("cp ~{sep = " " filters_bam_polyrad } .")
          system("cp ~{sep = " " filters_bam_supermassa } .")
          system("cp ~{sep = " " filters_bam_updog  } .")

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
            df_tot <- vector()
            for(j in 1:length(genotype_meth)){
              if(CountsFrom == "bam"){
                error.file <- paste0(SNPcall,"_error_bam_", genotype_meth[j],".txt")
              } else {
                error.file <- paste0(SNPcall,"_error_", genotype_meth[j],".txt")
              }
              error_df <- read.table(error.file, header = T, stringsAsFactors = F)
              error_df[,2] <- paste0("Chr10_",error_df[,2])
              colnames(error_df) <- c("ind", "mks", "gabGT", "methGT", "A", "AB", "BA", "B")
              df_meth <- merge(alleles, error_df)
              df_meth <- cbind(ErrorProb = genotype_meth[j], df_meth)
              df_tot <- rbind(df_tot, df_meth)
            }
            return(df_tot)
          }

          joint_maps <- function(CountsFrom = "vcf", genotype_meth = meth.geno, snpcall = "gatk"){
            map_df_tot <- coverage_df_tot<- vector()
            for(j in 1:length(genotype_meth)){
              if(CountsFrom == "vcf"){
                map.file <- paste0(snpcall, "_map_", genotype_meth[j],".txt")
              } else {
                map.file <- paste0(snpcall, "_map_bam_", genotype_meth[j],".txt")
              }

              map_df <- read.table(map.file, header = T)

              poscM <- (as.numeric(as.character(map_df[,2]))/1000000)*~{cmBymb}
              poscM.norm <- poscM-poscM[1]
              map_df <- cbind(map_df, poscM, poscM.norm)

              # Difference between distances real and estimated
              map_df <- cbind(map_df, diff= sqrt((map_df["poscM.norm"] - map_df["rf"])^2)[,1])

              coverage_df <- map_df[,2][length(map_df[,2])]*100/tot_mks[,2][length(tot_mks[,2])]
              map_df <- cbind(seed, depth, CountsFrom, ErrorProb= genotype_meth[j],
                              SNPcall=snpcall, map_df)
              coverage_df <- cbind(seed, depth, CountsFrom, ErrorProb= genotype_meth[j],
                                  SNPcall=snpcall, coverage = coverage_df)

              coverage_df_tot <- rbind(coverage_df_tot, coverage_df)
              map_df_tot <- rbind(map_df_tot, map_df)
            }
            return(list(map_df_tot, coverage_df_tot))
          }


            ########################################################################################
            # Table1: GenoCall; mks; ind; SNPcall; CountsFrom; alt; ref; gabGT; methGT; A; AB; BA; B
            ########################################################################################
            df_tot_all <- times <- maps_tot <- coverage_tot <- filters_tot <- vector()
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

            df_tot_bam <- joint_depths(alt.depth.bam,ref.depth.bam, "bam", method[i], depth, seed, meth.bam)

            ## Depths by softwares
            df_tot <- joint_depths(alt.depth,ref.depth, "vcf", method[i], depth, seed, meth.geno)

            df_tot_all <- rbind(df_tot_all, df_tot, df_tot_bam)

            ########################################################
            # Table2: seed; CountsFrom; ErrorProb; SNPcall; MK; rf; phases; real_phases
            # Table3: seed; CountsFrom; ErrorProb; SNPcall; coverage
            ########################################################

            map_df_tot <- joint_maps("vcf", meth.geno, method[i])
            map_bam_tot <- joint_maps("bam", meth.bam, method[i])
            map_df_gusmap <- joint_maps("vcf", "gusmap", method[i])
            map_bam_gusmap <- joint_maps("bam", "gusmap", method[i])

            maps_tot <- rbind(maps_tot, map_df_tot[[1]], map_bam_tot[[1]],
                              map_df_gusmap[[1]], map_bam_gusmap[[1]])
            maps_tot <- as.data.frame(maps_tot)
            coverage_tot <- rbind(coverage_tot, map_df_tot[[2]],
                                  map_bam_tot[[2]], map_df_gusmap[[2]],
                                  map_bam_gusmap[[2]])

            ##########################################################################
            # Table4: CountsFrom; seed; SNPcall; GenoCall; n_mks; distorted; redundant
            ##########################################################################

            filters_vcf_tot <- filters_bam_tot <- vector()
            for(j in 1:length(meth.filt)){
              filters_vcf <- read.table(paste0(method[i], "_filters_", meth.filt[j], ".txt"), header = T)
              filters_vcf <- cbind(CountsFrom = "vcf", seed= seed, depth, SNPcall = method[i], GenoCall = meth.filt[j],
                                  filters_vcf)
              filters_vcf_tot <- rbind(filters_vcf_tot, filters_vcf)
            }

            for(j in 1:length(meth.bam)){
              filters_bam <- read.table(paste0(method[i], "_filters_bam_", meth.bam[j], ".txt"), header = T)
              filters_bam <- cbind(CountsFrom = "bam", seed= seed, depth, SNPcall = method[i], GenoCall = meth.bam[j],
                                  filters_bam)
              filters_bam_tot <- rbind(filters_bam_tot, filters_bam)
            }

            filters_tot <- rbind(filters_tot, filters_bam_tot, filters_vcf_tot)

            ###########################################################################
            # Table6: CountsFrom; seed; SNPcall; GenoCall
            ###########################################################################

            CountsFrom <- vector()
            times_temp <- read.table(paste0(method[i], "_times.txt"), stringsAsFactors = F)
            temp <- strsplit(times_temp[,1], "_")
            temp <- do.call(rbind, temp)
            CountsFrom[which(temp[,3] == "bam")] <- "bam"
            CountsFrom[which(temp[,3] != "bam")] <- "vcf"
            SNPcall <- temp[,1]
            temp[which(temp[,3] == "bam"),3] <- temp[which(temp[,3] == "bam"),4]
            Genocall <- temp[,3]
            times <- rbind(times, data.frame(CountsFrom, seed, depth, SNPcall, Genocall, times = times_temp[,2]))
          }

          saveRDS(df_tot_all, file = "data1_depths_geno_prob.rds")
          saveRDS(maps_tot, file = "data2_maps.rds")
          saveRDS(coverage_tot, file = "data3_coverage.rds")
          saveRDS(filters_tot, file = "data4_filters.rds")
          saveRDS(times, file= "data6_times.rds")

      RSCRIPT
  >>>

  runtime{
    docker:"taniguti/onemap"
  }

  output{
    File data1_depths_geno_prob = "data1_depths_geno_prob.rds"
    File data2_maps = "data2_maps.rds"
    File data3_coverage = "data3_coverage.rds"
    File data4_filters = "data4_filters.rds"
    File data6_times   = "data6_times.rds"
  }
}
