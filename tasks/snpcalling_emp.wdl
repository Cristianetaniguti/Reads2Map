version 1.0

import "../structs/snpcalling_empS.wdl"
import "../structs/reference_struct.wdl"

import "create_alignment_from_families_files.wdl" as fam
import "gatk_genotyping.wdl" as gatk
import "freebayes_genotyping.wdl" as freebayes


workflow SNPCalling {

  input {
    Samples_info samples_info
    Reference references
    SplitVCF splitvcf
  }

  call fam.CreateAlignmentFromFamilies {
    input:
      families_info=samples_info.samples_info,
      references=references
  }

  call gatk.GatkGenotyping {
    input:
      alignments=CreateAlignmentFromFamilies.alignments,
      references=references,
      program="gatk",
      splitvcf = splitvcf,
      sampleNames = CreateAlignmentFromFamilies.names
  }

  call freebayes.FreebayesGenotyping {
    input:
      alignments=CreateAlignmentFromFamilies.alignments,
      bam=CreateAlignmentFromFamilies.bam,
      bai=CreateAlignmentFromFamilies.bai,
      references=references,
      program="freebayes",
      splitvcf = splitvcf,
      sampleNames = CreateAlignmentFromFamilies.names
  }

  output {
    File gatk_vcf_bi = GatkGenotyping.vcf_bi
    File gatk_vcf_multi = GatkGenotyping.vcf_multi
    File gatk_vcf_bi_bam_count = GatkGenotyping.vcf_bi_bam_counts
    File freebayes_vcf_bi = FreebayesGenotyping.vcf_bi
    File freebayes_vcf_multi = FreebayesGenotyping.vcf_multi
    File freebayes_vcf_bi_bam_count = FreebayesGenotyping.vcf_bi_bam_counts
  }
}
