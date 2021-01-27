version 1.0

import "reference_struct.wdl"

workflow CollectAllelicCountsToVcf {
    input {
        String program
        Array[String] sample_names
        Array[File] bams
        Array[File] bams_index
        Reference references
        File vcf_biallelics_splitted
        File vcf_biallelics_tbi_splitted
    }

    call CollectAllelicCounts {
        input:
        program=program,
        bams=bams,
        bams_index = bams_index,
        ref=references.ref_fasta,
        ref_fai=references.ref_fasta_index,
        ref_dict=references.ref_dict,
        vcf=vcf_biallelics_splitted,
        tbi=vcf_biallelics_tbi_splitted
    }

    call BamCounts4Onemap {
        input:
        sampleName=sample_names,
        counts=CollectAllelicCounts.counts,
        method = program
    }

    call BamDepths2Vcf {
        input:
        vcf_file = vcf_biallelics_splitted,
        ref_bam = BamCounts4Onemap.ref_bam,
        alt_bam = BamCounts4Onemap.alt_bam,
        example_alleles = BamCounts4Onemap.ref_alt_alleles,
        program = program
    }

    output {
        File vcf_biallelics_bamcounts = BamDepths2Vcf.bam_vcf
        File alt_bam = BamCounts4Onemap.alt_bam
        File ref_bam = BamCounts4Onemap.ref_bam
    }
}


task CollectAllelicCounts {
  input {
    String program
    Array[File] bams
    Array[File] bams_index
    File ref
    File ref_fai
    File ref_dict
    File vcf
    File tbi
  }

  Int disk_size = ceil(size(bams, "GB") + size(ref, "GB") * 2)

  command <<<
    set -e

    java -jar /gatk/picard.jar VcfToIntervalList \
      I=~{vcf} \
      O=interval.list

    mkdir links
    for bam in ~{sep=" " bams}; do ln -s $bam links/; done
    for bai in ~{sep=" " bams_index}; do ln -s $bai links/; done
    for bam in links/*.bam; do
      name=$(basename -s ".sorted.bam" "$bam")
      /gatk/gatk CollectAllelicCounts \
        --input "$bam" \
        --reference ~{ref} \
        --intervals interval.list \
        --output "${name}_~{program}_counts.tsv"
    done

  >>>

  runtime {
    docker: "taniguti/gatk-picard"
    memory: "4 GB"
    cpu: 1
    preemptible: 3
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    Array[File] counts = glob("*_~{program}_counts.tsv")
  }
}

# This task convert the output from BamCounts to the depths input for onemap
task BamCounts4Onemap{
  input {
    Array[File] counts
    Array[String] sampleName
    String method
  }

  Int disk_size = ceil(size(counts, "GB") + 5)

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

  runtime {
    docker: "cristaniguti/onemap_workflows"
    memory: "2 GB"
    cpu: 1
    preemptible: 4
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File ref_bam = "~{method}_ref_depth_bam.txt"
    File alt_bam = "~{method}_alt_depth_bam.txt"
    File ref_alt_alleles = "~{method}_ref_alt_alleles.txt"
  }
}

task BamDepths2Vcf{
  input {
    File vcf_file
    File ref_bam
    File alt_bam
    File example_alleles
    String program
  }

  Int disk_size = ceil(size(vcf_file, "GB") + size(ref_bam, "GB") + size(alt_bam, "GB") + 5)

  command <<<
    R --vanilla --no-save <<RSCRIPT

      library(onemap)
      library(vcfR)
      library(doParallel)
      source("/opt/scripts/functions_simu.R")

       ## Depths from bam
       depths.alt <- read.table("~{alt_bam}", header = T)
       depths.ref <- read.table("~{ref_bam}", header = T)

       depths <- list("ref" = depths.ref, "alt"=depths.alt)

       if(tail(strsplit("~{vcf_file}", "[.]")[[1]],1) =="gz") {
          vcf.temp <- paste0("vcf.temp",".", sample(1000,1), ".vcf")
          system(paste0("zcat ", "~{vcf_file}", " > ", vcf.temp))
          vcf_file <- vcf.temp
       } else {
          vcf_file <- "~{vcf_file}"
       }

       allele_file <- paste0("~{example_alleles}")
       bam_vcf <- make_vcf(vcf_file, depths, allele_file, "~{program}_bam_vcf.vcf")

       bam_vcfR <- read.vcfR(bam_vcf)
       save(bam_vcfR, file="~{program}_bam_vcfR.RData")

    RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/onemap_workflows"
    preemptible: 3
    memory: "4 GB"
    cpu: 3
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File bam_vcf = "~{program}_bam_vcf.vcf"
    File bam_vcfR = "~{program}_bam_vcfR.RData"
  }
}

