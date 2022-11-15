version 1.0

import "../../structs/dna_seq_structs.wdl"

import "../../subworkflows/create_alignment_from_families_files.wdl" as fam
import "../../subworkflows/gatk_genotyping.wdl" as gatk
import "../../subworkflows/freebayes_genotyping.wdl" as freebayes


workflow SNPCalling {

  input {
    File samples_info
    ReferenceFasta references
    Int max_cores
    Int chunk_size
    String rm_dupli
    String P1
    String P2
    String gatk_mchap
  }

  call fam.CreateAlignmentFromFamilies {
    input:
      families_info=samples_info,
      references=references,
      max_cores = max_cores,
      rm_dupli = rm_dupli,
      chunk_size = chunk_size
  }

  call gatk.GatkGenotyping {
    input:
      bams=CreateAlignmentFromFamilies.bam,
      bais=CreateAlignmentFromFamilies.bai,
      references=references,
      chunk_size = chunk_size,
      ploidy = 2,
      program="gatk",
      max_cores = max_cores,
      merged_bams = CreateAlignmentFromFamilies.merged_bam,
      P1 = P1,
      P2 = P2,
      mchap = gatk_mchap
  }

  call freebayes.FreebayesGenotyping {
    input:
      bams=CreateAlignmentFromFamilies.bam,
      bais=CreateAlignmentFromFamilies.bai,
      references=references,
      program="freebayes",
      max_cores = max_cores
  }

  output {
    File? gatk_multi_vcf = GatkGenotyping.vcf_multi
    File? gatk_vcf_bam_multi = GatkGenotyping.vcf_multi_bamcounts
    File gatk_vcf = GatkGenotyping.vcf_norm
    File gatk_vcf_bam_count = GatkGenotyping.vcf_norm_bamcounts
    File gatk_vcfEval = GatkGenotyping.vcfEval
    File Plots = GatkGenotyping.Plots
    File freebayes_vcf = FreebayesGenotyping.vcf_norm
    File freebayes_vcf_bam_count = FreebayesGenotyping.vcf_norm_bamcounts
    File freebayes_vcfEval = FreebayesGenotyping.vcfEval
    File merged_bam = CreateAlignmentFromFamilies.merged_bam
  }
}
