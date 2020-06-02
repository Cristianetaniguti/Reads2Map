version 1.0

import "../structs/alignment_struct.wdl"
import "../structs/reads_simuS.wdl"
import "../structs/snpcalling_empS.wdl"
import "./utils.wdl" as utils


workflow FreebayesGenotyping {
  input {
    Array[Alignment] alignments
    Array[File] bam
    Array[File] bai
    ReferenceFasta references
    Optional_filt optional_filt
    String program
  }

  call RunFreebayes {
    input:
      reference=references.ref_fasta,
      reference_idx=references.ref_fasta_index,
      bam=bam,
      bai=bai
  }

  call utils.TabixVcf {
    input:
      variants=RunFreebayes.vcf
  }

  call utils.VcftoolsApplyFilters {
    input:
      vcf_in=TabixVcf.vcf,
      max_missing=0.75,
      min_alleles=2,
      max_alleles=2,
      maf=optional_filt.maf,
      program=program,
      min_meanDP = optional_filt.min_meanDP
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
    File reference
    File reference_idx
    Array[File] bam
    Array[File] bai
  }

  command <<<
   freebayes-parallel <(fasta_generate_regions.py ~{reference_idx} 100000) 20 \
   --genotype-qualities -f ~{reference}  ~{sep=" " bam} > "freebayes.vcf"
  >>>

  runtime {
    docker: "taniguti/freebayes"
    mem:"--nodes=1"
    time:"24:00:00"
    cpu:20
  }

  output {
    File vcf = "freebayes.vcf"
  }
}
