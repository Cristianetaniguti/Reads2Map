version 1.0

import "../structs/alignment_struct.wdl"
import "../structs/reads_simuS.wdl"
import "./utils.wdl" as utils


workflow FreebayesGenotyping {
  input {
    Alignment alignment
    ReferenceFasta references

  }

  call RunFreebayes {
    input:
      sample=alignment.sample,
      reference=references.ref_fasta,
      bam=alignment.bam,
      bai=alignment.bai
  }

  call utils.TabixVcf {
    input:
      sample=alignment.sample,
      variants=RunFreebayes.vcf
  }

  call utils.BamCounts {
    input:
      sample=alignment.sample,
      program="freebayes",
      bam=alignment.bam,
      bai=alignment.bai,
      ref=references.ref_fasta,
      ref_fai=references.ref_fasta_index,
      ref_dict=references.ref_dict,
      vcf=TabixVcf.vcf,  # Precisa ser o VCF consolidado
      tbi=TabixVcf.tbi
    }

  output {
    File vcf = TabixVcf.vcf
    File tbi = TabixVcf.tbi
    File counts = BamCounts.counts
  }
}


# Variant calling using freebayes
task RunFreebayes {

  input {
    String sample
    File reference
    File bam
    File bai
  }

  command <<<
    freebayes --genotype-qualities -f ~{reference} ~{bam} > ~{sample}.freebayes.vcf
  >>>

  runtime {
    docker: "taniguti/freebayes"
  }

  output {
    File vcf = "~{sample}.freebayes.vcf"
  }
}
