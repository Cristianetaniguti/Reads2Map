version 1.0

import "../structs/dna_seq_structs.wdl"
import "../tasks/utils.wdl" as utils
import "../tasks/freebayes.wdl"
import "../tasks/chunk_lists.wdl"

import "../subworkflows/norm_filt_vcf.wdl" as norm_filt


workflow FreebayesGenotyping {
  input {
    File merged_bam
    #Array[File] bams
    #Array[File] bais
    ReferenceFasta references
    String program
    Int max_cores
    Int ploidy
    File? vcf_simu
    Boolean replaceAD
    Int n_chrom
  }

  call chunk_lists.CreateChunksBamByChr {
    input:
      merged_bam = merged_bam,
      reference = references.ref_fasta,
      n_chrom = n_chrom
  }

  scatter (chunk in zip(CreateChunksBamByChr.bams_chunks, CreateChunksBamByChr.bais_chunks)) {

    call freebayes.RunFreebayes {
      input:
        reference = references.ref_fasta,
        reference_idx=references.ref_fasta_index,
        bam = chunk.left,
        bai = chunk.right,
        max_cores = max_cores,
        ploidy = ploidy
    }
  }

  call utils.mergeVCFs {
    input:
      haplo_vcf = RunFreebayes.vcf
  }

  call norm_filt.Normalization {
    input:
      vcf_in = mergeVCFs.merged_vcf,
      vcf_simu = vcf_simu,
      reference = references.ref_fasta,
      reference_idx = references.ref_fasta_index,
      reference_dict = references.ref_dict
  }

  Map[String, Array[File]] map_bams = {"bam": CreateChunksBamByChr.bams_chunks, "bai": CreateChunksBamByChr.bais_chunks}

  if(replaceAD){
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
  }

  output {
    File vcf_norm = Normalization.vcf_norm
    File? vcf_norm_bamcounts = ReplaceAD.bam_vcf
    File vcfEval = Normalization.vcfEval
  }
}
