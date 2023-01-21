version 1.0

import "../tasks/bcftools.wdl"
import "../tasks/gatk.wdl"

workflow Normalization {
  input {
    File vcf_in
    File? vcf_simu
    File reference
    File reference_idx
    File reference_dict
    String program
    String counts_source
    Int ploidy
  }

  call bcftools.BiallelicNormalization {
    input:
      vcf_file = vcf_in,
      reference = reference,
      reference_idx = reference_idx,
      ploidy = ploidy,
      software = program
  }

  call gatk.VariantEval {
    input:
      vcf_norm = BiallelicNormalization.vcf_norm,
      vcf_norm_tbi = BiallelicNormalization.vcf_norm_tbi,
      vcf_simu = vcf_simu,
      reference = reference,
      reference_idx = reference_idx,
      reference_dict = reference_dict,
      ploidy = ploidy
  }

  output {
    File vcf_norm = BiallelicNormalization.vcf_norm
    File vcf_norm_tbi = BiallelicNormalization.vcf_norm_tbi
    File vcfEval = VariantEval.vcfEval
    String software = "~{program}"
    String source = "~{counts_source}"
  }
}
