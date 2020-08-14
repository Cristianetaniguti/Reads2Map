version 1.0

import "./tasks/maps_emp.wdl" as maps


workflow EmpiricalReads {

  input {
    File gatk_vcf
    File freebayes_vcf
    File gatk_vcf_bi_bam_count
    File freebayes_vcf_bi_bam_count
    Dataset dataset
    String? filters
  }

  call maps.Maps {
    input:
      dataset = dataset,
      gatk_vcf = gatk_vcf,
      freebayes_vcf = freebayes_vcf,
      gatk_vcf_bam_counts = gatk_vcf_bi_bam_count,
      freebayes_vcf_bam_counts = freebayes_vcf_bi_bam_count,
      filters = filters
  }


  output {
    File EmpiricalReads_results = Maps.EmpiricalReads_results
  }
}
