version 1.0

import "../structs/dna_seq_structs.wdl"
import "../tasks/utils.wdl" as utils
import "../tasks/freebayes.wdl"

import "norm_filt_vcf.wdl" as norm_filt


workflow FreebayesGenotyping {
  input {
    Array[File] bams
    Array[File] bais
    ReferenceFasta references
    String program
    Int max_cores
    File? vcf_simu
  }

  call freebayes.RunFreebayes {
    input:
      reference=references.ref_fasta,
      reference_idx=references.ref_fasta_index,
      bam=bams,
      bai=bais,
      max_cores = max_cores
  }

  call norm_filt.Normalization {
    input:
      vcf_in = RunFreebayes.vcf,
      vcf_simu = vcf_simu,
      reference = references.ref_fasta,
      reference_idx = references.ref_fasta_index,
      reference_dict = references.ref_dict
  }

  Map[String, Array[File]] map_bams = {"bam": bams, "bai": bais}

  call utils.ReplaceAD {
    input:
      ref_fasta = references.ref_fasta,
      ref_index = references.ref_fasta_index,
      bams = map_bams["bam"],
      bais = map_bams["bai"],
      vcf = Normalization.vcf_norm,
      tbi = Normalization.vcf_norm_tbi,
      program = program
  }

  output {
    File vcf_norm = Normalization.vcf_norm
    File vcf_norm_bamcounts = ReplaceAD.bam_vcf
    File vcfEval = Normalization.vcfEval
  }
}
