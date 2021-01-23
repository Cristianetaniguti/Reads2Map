version 1.0

import "../structs/alignment_struct.wdl"

task TabixVcf {
  input {
    File variants
  }

  command <<<
    bgzip -c ~{variants} > freebayes.vcf.gz
    tabix -p vcf freebayes.vcf.gz

  >>>

  runtime {
    docker: "taniguti/gatk-picard"
    time:"24:00:00"
    mem:"--nodes=1"
    cpu:1
  }

  output {
    File vcf = "freebayes.vcf.gz"
    File tbi = "freebayes.vcf.gz.tbi"
  }
}


task BamCounts {
  input {
    String program
    Array[File] bam
    File ref
    File ref_fai
    File ref_dict
    File vcf
    File tbi
  }

  command <<<
    set -e

    java -jar /gatk/picard.jar VcfToIntervalList \
      I=~{vcf} \
      O=interval.list

    for file in ~{sep= " " bam}; do

      samtools index $file
      sample=`basename -s .sorted.bam $file`
      echo $sample

      /gatk/gatk CollectAllelicCounts \
        --input $file \
        --reference ~{ref} \
        --intervals interval.list \
        --output "${sample}_~{program}_counts.tsv"

    done

  >>>

  runtime{
    docker:"taniguti/gatk-picard"
    mem:"--nodes=1"
    cpu:1
    time:"48:00:00"
  }

  output{
    Array[File] counts = glob("*_~{program}_counts.tsv")
  }
}


task VcftoolsMerge {

  input {
    String prefix
    Array[File] vcfs
    Array[File] tbis
  }

  command <<<
    echo "~{sep=' ' tbis}"
    vcf-merge ~{sep=" "  vcfs} > ~{prefix}.variants.vcf
    bgzip ~{prefix}.variants.vcf
    tabix -p vcf ~{prefix}.variants.vcf.gz

  >>>
  runtime {
    docker: "taniguti/vcftools"
    mem:"--nodes=1"
    cpu:1
    time:"24:00:00"
  }

  output {
    File vcf = "~{prefix}.variants.vcf.gz"
    File tbi = "~{prefix}.variants.vcf.gz.tbi"
  }
}

task BcftoolsMerge {

  input {
    String prefix
    Array[File] vcfs
    Array[File] tbis
  }

  command <<<
    echo "~{sep=' ' tbis}"
    bcftools merge ~{sep=" "  vcfs} > ~{prefix}.variants.vcf

  >>>
  runtime {
    docker: "biocontainers/bcftools:1.3.1"
    mem:"--mem-per-cpu=24042"
    cpu:1
    time:"48:00:00"
  }

  output {
    File vcf = "~{prefix}.variants.vcf"
  }
}


task VcftoolsApplyFilters {

  input {
    File vcf_in
    Float max_missing
    Int min_alleles
    Int max_alleles
    Float? maf
    String program
    Int? min_meanDP
    String? chromosome
  }

  command <<<
    vcftools --gzvcf ~{vcf_in} --max-missing ~{max_missing} --min-alleles ~{min_alleles} --max-alleles ~{max_alleles} ~{"--maf " +  maf} ~{"--min-meanDP " +  min_meanDP} ~{"--chr " +  chromosome} --recode --out ~{program}

    bgzip ~{program}.recode.vcf
    tabix -p vcf ~{program}.recode.vcf.gz

  >>>
  runtime {
    docker: "taniguti/vcftools"
    mem:"--nodes=1"
    cpu:1
    time:"24:00:00"
  }

  output {
    File vcf = "~{program}.recode.vcf.gz"
    File tbi = "~{program}.recode.vcf.gz.tbi"
  }
}


task CalculateVcfMetrics {

  input {
    File freebayesVCF
    File gatkVCF
    File ref_alt_alleles
    Int seed
    Int depth
  }

  command <<<

        R --vanilla --no-save <<RSCRIPT

        # If variants are simulated by pirs,
        # the metrics will be related to the total simulated not only the captured by RAD
        library(vcfR)
        freebayes <- read.vcfR("~{freebayesVCF}")
        gatk <- read.vcfR("~{gatkVCF}")

        snps <- read.table("~{ref_alt_alleles}", stringsAsFactors = F)
        simulated.pos <- snps[,2]
        simulated.ref <- snps[,3]
        simulated.alt <- snps[,4]

        methods <- c("freebayes", "gatk")
        results_tot <- vector()

        for(i in methods){
          # counting corrected identified markers
          pos <- as.numeric(as.character(get(i)@fix[,2]))
          chr <- get(i)@fix[,1]
          site_list <- data.frame(chr, pos, pos)
          # Export for next step
          write.table(site_list, file= paste0(i,"_site_list.txt"), quote=F, row.names=F, sep="\t", col.names=F)

          ref <- get(i)@fix[,4]
          alt <- get(i)@fix[,5]

          nmk.filt <- length(simulated.pos)
          nmk.id <- length(pos)

          ok <- sum(simulated.pos %in% pos)
          falso.positivo <- sum(!(pos %in% simulated.pos))
          ref.ok <- sum(simulated.ref==ref[pos %in% simulated.pos])
          alt.ok <- sum(simulated.alt==alt[pos %in% simulated.pos])

          result <- data.frame(depth = ~{depth}, seed = ~{seed}, SNPcall = i,mks_tot = nmk.filt, mks_ide = nmk.id, ok, fake=falso.positivo, ref.ok, alt.ok)
          results_tot <- rbind(results_tot, result)

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
    docker: "cristaniguti/onemap_workflows"
    mem:"--nodes=1"
    cpu:1
    time:"24:00:00"
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


# This task convert the output from BamCounts to the depths input for onemap
task BamCounts4Onemap{
  input{
    Array[File] counts
    Array[String] sampleName
    String method
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(R.utils)
      system("cp ~{sep=" "  counts} .")
      names <- c("~{sep=" , "  sampleName}")
      names <- unlist(strsplit(names, split = " , "))

      system(paste("grep -n 'CONTIG'", paste0(names[1],"_", "~{method}","_counts.tsv"), "> idx"))
      idx <- read.table("idx", stringsAsFactors = F)
      idx <- strsplit(idx[,1], ":")[[1]][1]
      idx <- as.numeric(idx) -1

      file.counts <- read.table(paste0(names[1],"_", "~{method}","_counts.tsv"), skip = idx, header=T, stringsAsFactors = F)

      ref_depth_matrix2 <- alt_depth_matrix2  <- matrix(NA, nrow = dim(file.counts)[1], ncol = length(names))

      for(j in 1:length(names)){
        ## From picard tool

        file.counts <- read.table(paste0(names[j],"_", "~{method}","_counts.tsv"), skip = idx, header=T, stringsAsFactors = F)

        ref_depth_matrix2[,j] <- file.counts[,3]
        alt_depth_matrix2[,j] <- file.counts[,4]

        if (j == 1){
          ref_allele <- file.counts[,5]
          alt_allele <- file.counts[,6]
        } else {
          idx.ref <- which(ref_allele == "N")
          idx.alt <- which(alt_allele == "N")
          if (length(idx.ref)!=0){
            ref_allele[idx.ref] <- file.counts[idx.ref,5]
          }
          if (length(idx.alt)!=0){
            alt_allele[idx.alt] <- file.counts[idx.alt,6]
          }
        }

        rownames(ref_depth_matrix2) <- rownames(alt_depth_matrix2) <- paste0(file.counts[,1],"_", file.counts[,2])
        colnames(ref_depth_matrix2) <- colnames(alt_depth_matrix2) <- names

        alleles <- data.frame(file.counts[,1],file.counts[,2], ref_allele, alt_allele)
        write.table(alleles, file = paste0("~{method}","_ref_alt_alleles.txt"), col.names = F, row.names = F)

        write.table(ref_depth_matrix2, file = paste0("~{method}","_ref_depth_bam.txt"), quote=F, row.names=T, sep="\t", col.names=T)
        write.table(alt_depth_matrix2, file = paste0("~{method}","_alt_depth_bam.txt"), quote=F, row.names=T, sep="\t", col.names=T)
      }

    RSCRIPT

  >>>

  runtime{
    docker:"cristaniguti/onemap_workflows"
    mem:"--nodes=1"
    cpu:1
    time:"48:00:00"
  }

  output{
    File ref_bam      = "~{method}_ref_depth_bam.txt"
    File alt_bam      = "~{method}_alt_depth_bam.txt"
    File ref_alt_alleles    = "~{method}_ref_alt_alleles.txt"
  }
}

task ApplyRandomFilters{
  input{
    File gatk_vcf
    File freebayes_vcf
    File gatk_vcf_bam_counts
    File freebayes_vcf_bam_counts
    String? filters
    String? chromosome
  }

  command <<<
    vcftools --gzvcf ~{gatk_vcf} ~{filters} ~{"--chr " + chromosome} --recode --stdout > gatk_vcf_filt.vcf

    vcftools --gzvcf ~{freebayes_vcf} ~{filters} ~{"--chr " + chromosome} --recode --stdout > freebayes_vcf_filt.vcf

    vcftools --gzvcf ~{gatk_vcf_bam_counts} ~{filters} ~{"--chr " + chromosome} --recode --stdout > gatk_vcf_bam_counts_filt.vcf

    vcftools --gzvcf ~{freebayes_vcf_bam_counts} ~{filters} ~{"--chr " + chromosome} --recode --stdout > freebayes_vcf_bam_counts_filt.vcf
  >>>

  runtime{
    docker:"taniguti/vcftools"
    mem:"--nodes=1"
    cpu:1
    time:"24:00:00"
  }

  output{
    File gatk_vcf_filt = "gatk_vcf_filt.vcf"
    File freebayes_vcf_filt = "freebayes_vcf_filt.vcf"
    File gatk_vcf_bam_counts_filt = "gatk_vcf_bam_counts_filt.vcf"
    File freebayes_vcf_bam_counts_filt = "freebayes_vcf_bam_counts_filt.vcf"
  }
}

task Gambis {
  input{
    File gatk_vcf_bam_counts
    File freebayes_vcf_bam_counts
    String method
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
      if("~{method}" == "gatk") choosed_bam = "~{gatk_vcf_bam_counts}" else choosed_bam = "~{freebayes_vcf_bam_counts}"

      system(paste("cp", choosed_bam, "choosed_bam.vcf"))

    RSCRIPT
  >>>

  output{
    File choosed_bam = "choosed_bam.vcf"
  }
}

task FiltChr {
  input {
    File vcf_bam
    String chromosome
  }

  command <<<
    vcftools --gzvcf ~{vcf_bam} --chr ~{chromosome} --recode --stdout > bam_chr.vcf
  >>>

  runtime {
    docker:"taniguti/vcftools"
    mem:"--nodes=1"
    cpu:1
    time:"24:00:00"
  }

  output {
    File bam_chr = "bam_chr.vcf"
  }
}
