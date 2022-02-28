version 1.0

import "structs/struct_reference.wdl"

import "tasks/create_alignment_from_families_files.wdl" as fam
import "tasks/gatk_genotyping.wdl" as gatk


workflow SNPCalling_gatk {

  input {
    File samples_info
    Reference references
    Int max_cores
    Int chunk_size
    Int ploidy
    String rm_dupli
    String P1
    String P2
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
      ploidy = ploidy,
      program="gatk",
      P1 = P1,
      P2 = P2
  }

  output {
    File gatk_vcf = GatkGenotyping.vcf_norm
    File gatk_vcf_bam_count = GatkGenotyping.vcf_norm_bamcounts
    File gatk_vcfEval = GatkGenotyping.vcfEval
    File Plots = GatkGenotyping.Plots
    File merged_bam = CreateAlignmentFromFamilies.merged_bam
  }
}
