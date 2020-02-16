version 1.0

import "../structs/alignment_struct.wdl"
import "../structs/reads_simuS.wdl"
import "./utils.wdl" as utils


workflow FreebayesGenotyping {
  input {
    Array[Alignment] alignments
    ReferenceFasta references
    String program

  }

  scatter (alignment in alignments) {
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
  }

  call utils.BcftoolsMerge {
    input:
      prefix=program,
      vcfs=TabixVcf.vcf,
      tbis=TabixVcf.tbi
  }

  call utils.VcftoolsApplyFilters {
    input:
      vcf_in=BcftoolsMerge.vcf,
      max_missing=0.75,
      min_alleles=2,
      max_alleles=2,
      maf=0.05,
      program=program
  }


  scatter (alignment in alignments) {
    call utils.BamCounts {
      input:
        sample=alignment.sample,
        program=program,
        bam=alignment.bam,
        bai=alignment.bai,
        ref=references.ref_fasta,
        ref_fai=references.ref_fasta_index,
        ref_dict=references.ref_dict,
        vcf=VcftoolsApplyFilters.vcf,
        tbi=VcftoolsApplyFilters.tbi
    }
  }

  output {
    File vcf = VcftoolsApplyFilters.vcf
    File tbi = VcftoolsApplyFilters.tbi
    Array[File] counts = BamCounts.counts
  }
}


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
