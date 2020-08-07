version 1.0

import "../structs/alignment_struct.wdl"
import "../structs/reads_simuS.wdl"
import "../structs/snpcalling_empS.wdl"
import "./utils.wdl" as utils
import "split_filt_vcf.wdl" as norm_filt


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

  call norm_filt.SplitFiltVCF{
    input:
      vcf_in=RunFreebayes.vcf,
      max_missing=0.75,
      maf=optional_filt.maf,
      program=program,
      min_meanDP = optional_filt.min_meanDP,
      chromosome = optional_filt.chromosome,
      reference = references.ref_fasta,
      reference_idx = references.ref_fasta_index,
      parent1 = optional_filt.parent1,
      parent2 = optional_filt.parent2
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
        vcf=SplitFiltVCF.vcf_bi_chr_norm,
        tbi=SplitFiltVCF.vcf_bi_chr_norm_tbi

    }
  }

  output {
    File vcf = SplitFiltVCF.vcf_bi_chr_norm
    File tbi = SplitFiltVCF.vcf_bi_chr_norm_tbi
    Array[File] counts = BamCounts.counts
    File vcf_bi_tot = SplitFiltVCF.vcf_bi_norm
    File vcf_multi_tot = SplitFiltVCF.vcf_multi_norm
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
    time:"72:00:00"
    cpu:20
  }

  output {
    File vcf = "freebayes.vcf"
  }
}
