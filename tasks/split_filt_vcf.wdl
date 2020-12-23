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
    File vcf_bi = SplitFilters.vcf_bi
    File vcf_bi_tbi = SplitFilters.vcf_bi_tbi
    File vcf_multi = SplitFilters.vcf_multi
  }
}


task BiallelicNormalization {
  input {
    File vcf_file
    File reference
    File reference_idx
  }

  command <<<
    bcftools norm ~{vcf_file} --rm-dup all -Ov --check-ref w -f ~{reference} > vcf_norm.vcf
  >>>

  runtime {
    docker: "lifebitai/bcftools"
    mem:"--nodes=1"
    time:"72:00:00"
    cpu:1
  }

  output {
    File vcf_norm = "vcf_norm.vcf"
  }
}


task SplitFilters{
  input{
    File vcf_in
    String program
    String parent1
    String parent2
  }

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

  >>>

  runtime {
    docker: "cristaniguti/vcftools"
    mem:"--nodes=1"
    time:"72:00:00"
    cpu:1
  }

  output {
    File vcf_multi = "~{program}_multi.recode.vcf.gz"
    File vcf_bi = "~{program}_bi.recode.vcf.gz"
    File vcf_bi_tbi = "~{program}_bi.recode.vcf.gz.tbi"
  }
}


