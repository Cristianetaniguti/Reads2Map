version 1.0

workflow SplitFiltVCF{
  input {
    File vcf_in
    String program
    File reference
    File reference_idx
    String parent1
    String parent2
  }

  call BiallelicNormalization {
    input:
      vcf_file = vcf_in,
      reference = reference,
      reference_idx = reference_idx
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
    memory: "2 GB"
    preemptible: 3
    cpu: 1
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File vcf_norm = "vcf_norm.vcf"
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
    vcftools --gzvcf ~{vcf_in}  --min-alleles 3 --recode --out ~{program}_multi1

    vcftools --gzvcf ~{vcf_in}  --min-alleles 2 --max-alleles 2 --recode --out ~{program}_bi1

    Rscript /opt/scripts/split.R ~{program}_bi1.recode.vcf ~{parent1} ~{parent2} position_multi2.txt

    vcftools --vcf ~{program}_bi1.recode.vcf --positions position_multi2.txt --recode --out ~{program}_multi2

    vcftools --vcf ~{program}_bi1.recode.vcf --exclude-positions position_multi2.txt --recode --out ~{program}_bi

    vcf-concat ~{program}_multi1.recode.vcf ~{program}_multi2.recode.vcf  > ~{program}_multi.recode.vcf

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


