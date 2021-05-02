version 1.0

import "structs/snpcalling_empS.wdl"
import "structs/reference_struct.wdl"

import "tasks/create_alignment_from_families_files.wdl" as fam
import "tasks/gatk_genotyping.wdl" as gatk
import "tasks/freebayes_genotyping.wdl" as freebayes


workflow SNPCalling {

  input {
    Samples_info samples_info
    Reference references
    SplitVCF splitvcf
    Int max_cores
  }

  call fam.CreateAlignmentFromFamilies {
    input:
      families_info=samples_info.samples_info,
      references=references,
      max_cores = max_cores
  }

  call gatk.GatkGenotyping {
    input:
      bams=CreateAlignmentFromFamilies.bam,
      bais=CreateAlignmentFromFamilies.bai,
      references=references,
      program="gatk",
      parent1 = splitvcf.parent1,
      parent2 = splitvcf.parent2
  }

  call freebayes.FreebayesGenotyping {
    input:
      bams=CreateAlignmentFromFamilies.bam,
      bais=CreateAlignmentFromFamilies.bai,
      references=references,
      program="freebayes",
      parent1 = splitvcf.parent1,
      parent2 = splitvcf.parent2,
      max_cores = max_cores
  }

  output {
    File gatk_vcf_bi = GatkGenotyping.vcf_biallelics
    File gatk_vcf_multi = GatkGenotyping.vcf_multiallelics
    File gatk_vcf_bi_bam_count = GatkGenotyping.vcf_biallelics_bamcounts
    File gatk_vcfEval = GatkGenotyping.vcfEval
    File Plots = GatkGenotyping.Plots
    File freebayes_vcf_bi = FreebayesGenotyping.vcf_biallelics
    File freebayes_vcf_multi = FreebayesGenotyping.vcf_multiallelics
    File freebayes_vcf_bi_bam_count = FreebayesGenotyping.vcf_biallelics_bamcounts
    File freebayes_vcfEval = FreebayesGenotyping.vcfEval
  }
}
