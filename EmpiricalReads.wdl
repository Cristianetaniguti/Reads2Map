version 1.0

import "structs/snpcalling_empS.wdl"
import "structs/maps_empS.wdl"

import "tasks/snpcalling_emp.wdl" as snpcalling
import "tasks/maps_emp.wdl" as maps

workflow EmpiricalReads {

  input {
    Samples_info samples_info
    References references
    Dataset dataset
  }

  call snpcalling.SNPCalling{
    input:
      samples_info = samples_info,
      references = references
  }

  call maps.Maps{
    input:
      dataset = dataset,
      gatk_vcf = SNPCalling.gatk_vcf,
      freebayes_vcf = SNPCalling.freebayes_vcf,
      freebayes_ref_bam = SNPCalling.freebayes_ref_bam,
      freebayes_alt_bam = SNPCalling.freebayes_alt_bam,
      gatk_ref_bam = SNPCalling.gatk_ref_bam,
      gatk_alt_bam = SNPCalling.gatk_alt_bam,
      gatk_example_alleles = SNPCalling.gatk_example_alleles,
      freebayes_example_alleles = SNPCalling.freebayes_example_alleles
  }
  
  output{
    File EmpiricalReads_results = Maps.EmpiricalReads_results
    File gatk_vcf = SNPCalling.gatk_vcf_tot
    File freebayes_vcf = SNPCalling.freebayes_vcf_tot
  }
}
