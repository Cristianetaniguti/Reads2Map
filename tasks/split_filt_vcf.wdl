version 1.0

workflow SplitFiltVCF{
  input {
    File vcf_in
    File? vcf_simu
    String program
    File reference
    File reference_idx
    File reference_dict
    String parent1
    String parent2
  }

  call BiallelicNormalization {
    input:
      vcf_file = vcf_in,
      reference = reference,
      reference_idx = reference_idx
  }

  call VariantEval {
    input:
      vcf_norm = BiallelicNormalization.vcf_norm,
      vcf_simu = vcf_simu,
      reference = reference,
      reference_idx = reference_idx,
      reference_dict = reference_dict
  }

  call SplitFilters {
    input:
      vcf_in = BiallelicNormalization.vcf_norm,
      program = program,
      parent1 = parent1,
      parent2 = parent2
  }

  output {
    File vcf_biallelics = SplitFilters.vcf_biallelics
    File vcf_biallelics_tbi = SplitFilters.vcf_biallelics_tbi
    File vcf_multiallelics = SplitFilters.vcf_multiallelics
    File vcfEval = VariantEval.vcfEval
  }
}

task BiallelicNormalization {
  input {
    File vcf_file
    File reference
    File reference_idx
  }

  Int disk_size = ceil(size(vcf_file, "GB") + size(reference, "GB") + 2)

  command <<<
    bcftools norm ~{vcf_file} --rm-dup all -Ov --check-ref w -f ~{reference} > vcf_norm.vcf
  >>>

  runtime {
    docker: "lifebitai/bcftools"
    memory: "1 GB"
    preemptible: 3
    cpu: 1
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File vcf_norm = "vcf_norm.vcf"
  }
}

task VariantEval {
  input {
    File vcf_norm
    File? vcf_simu
    File reference
    File reference_idx
    File reference_dict
  }

  Int disk_size = ceil(size(vcf_norm, "GB") + size(reference, "GB") + size(vcf_simu, "GB") + 2)

  command <<<
    java -jar /usr/GenomeAnalysisTK.jar -T VariantEval -R ~{reference} -eval ~{vcf_norm} ~{"-D " + vcf_simu} -EV ValidationReport -EV CountVariants -o vcfEval.txt
  >>>

  runtime {
    docker: "broadinstitute/gatk3:3.8-1"
    memory: "1 GB"
    preemptible: 3
    cpu: 1
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File vcfEval = "vcfEval.txt"
  }

}

task SplitFilters {
  input {
    File vcf_in
    String program
    String parent1
    String parent2
  }

  Int disk_size = ceil(size(vcf_in, "GB") + 2)

  command <<<
    # TODO: Change for bcftools
    vcftools --gzvcf ~{vcf_in}  --min-alleles 3 --recode --out ~{program}_multi1

    vcftools --gzvcf ~{vcf_in}  --min-alleles 2 --max-alleles 2 --recode --out ~{program}_bi1

    Rscript /opt/scripts/split.R ~{program}_bi1.recode.vcf ~{parent1} ~{parent2} position_multi2.txt

    vcftools --vcf ~{program}_bi1.recode.vcf --positions position_multi2.txt --recode --out ~{program}_multi2

    vcftools --vcf ~{program}_bi1.recode.vcf --exclude-positions position_multi2.txt --recode --out ~{program}_bi

    vcf-concat ~{program}_multi1.recode.vcf ~{program}_multi2.recode.vcf  > ~{program}_multi.recode.vcf

    vcftools --vcf ~{program}_multi.recode.vcf --maf 0.05 --out ~{program}_multi2 --recode

    mv ~{program}_multi2.recode.vcf ~{program}_multi.recode.vcf

   # sort
    grep "^#" ~{program}_multi.recode.vcf > output.vcf
    grep -v "^#" ~{program}_multi.recode.vcf| sort -k1,1V -k2,2g >> output.vcf

    mv output.vcf ~{program}_multi.recode.vcf

    bgzip ~{program}_multi.recode.vcf    

    bgzip ~{program}_bi.recode.vcf
    tabix -p vcf ~{program}_bi.recode.vcf.gz
    tabix -p vcf ~{program}_multi.recode.vcf.gz

  >>>

  runtime {
    docker: "cristaniguti/vcftools"
    memory: "2 GB"
    preemptible: 3
    cpu: 1
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File vcf_multiallelics = "~{program}_multi.recode.vcf.gz"
    File vcf_multiallelics_tbi = "~{program}_multi.recode.vcf.gz.tbi"
    File vcf_biallelics = "~{program}_bi.recode.vcf.gz"
    File vcf_biallelics_tbi = "~{program}_bi.recode.vcf.gz.tbi"
  }
}


