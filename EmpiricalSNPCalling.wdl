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
    Int max_cores
    Int chunk_size
    String rm_dupli
  }

  call fam.CreateAlignmentFromFamilies {
    input:
      families_info=samples_info.samples_info,
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
      program="gatk"
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
