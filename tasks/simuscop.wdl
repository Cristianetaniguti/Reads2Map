version 1.0

import "../structs/dna_seq_structs.wdl"

task SimuscopProfile {
  input {
    String library_type
    File?  emp_bam
    File   vcf
    ReferenceFasta  references
  }

  Int disk_size = ceil(size(vcf, "GiB") * 2 + size(emp_bam, "GiB"))
  Int memory_size = 5000

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
    docker: "cristaniguti/reads2map:0.0.4"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "SimuscopProfile"
    mem:"~{memory_size}M"
    time:"05:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Run [seqToProfile](https://github.com/qasimyu/simuscop) to generate data set profile."
  }

  output {
    File profile = "sample.profile"
  }
}

task SimuscopSimulation {
 input {
    String library_type
    String sampleName
    Int depth
    File? emp_bam
    File vcf
    ReferenceFasta references
    String chrom
    File profile
  }

  Int disk_size = ceil(size(emp_bam, "GiB") + size(vcf, "GiB") + size(references.ref_fasta, "GiB") * depth)
  Int memory_size = 10000

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
    docker: "cristaniguti/reads2map:0.0.4"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "SimuscopSimulation"
    mem:"~{memory_size}M"
    time:"10:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Run [simuReads](https://github.com/qasimyu/simuscop) to simulated exome or WGS sequencing reads."
  }

  output {
    File fastq_seq = "~{sampleName}.fq"
  }
}