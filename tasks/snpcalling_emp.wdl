version 1.0

import "../structs/snpcalling_empS.wdl"

import "create_alignment_from_families_files.wdl" as fam
import "gatk_genotyping.wdl" as gatk
import "freebayes_genotyping.wdl" as freebayes
import "utils.wdl" as utils
import "utilsR.wdl" as utilsR

workflow SNPCalling{

  input {
    Samples_info samples_info
    References references
    Optional_filt optional_filt
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
      optional_filt = optional_filt
  }

  call freebayes.FreebayesGenotyping {
    input:
      alignments=CreateAlignmentFromFamilies.alignments,
      bam=CreateAlignmentFromFamilies.bam,
      bai=CreateAlignmentFromFamilies.bai,
      references=references,
      program="freebayes",
      optional_filt = optional_filt
  }

  call utils.BamCounts4Onemap {
    input:
      sampleName=CreateAlignmentFromFamilies.names,
      freebayes_counts=FreebayesGenotyping.counts,
      gatk_counts=GatkGenotyping.counts
  }

  output{
    File gatk_vcf = GatkGenotyping.vcf
    File gatk_vcf_bi_tot = GatkGenotyping.vcf_bi_tot
    File gatk_vcf_multi_tot = GatkGenotyping.vcf_multi_tot
    File freebayes_vcf = FreebayesGenotyping.vcf
    File freebayes_vcf_bi_tot = FreebayesGenotyping.vcf_bi_tot
    File freebayes_vcf_multi_tot = FreebayesGenotyping.vcf_multi_tot
    File freebayes_ref_bam = BamCounts4Onemap.freebayes_ref_bam
    File freebayes_alt_bam = BamCounts4Onemap.freebayes_alt_bam
    File gatk_ref_bam = BamCounts4Onemap.gatk_ref_bam
    File gatk_alt_bam = BamCounts4Onemap.gatk_alt_bam
    File gatk_example_alleles = BamCounts4Onemap.gatk_example_alleles
    File freebayes_example_alleles = BamCounts4Onemap.freebayes_example_alleles
  }
}