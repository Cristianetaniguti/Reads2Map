version 1.0

workflow Normalization{
  input {
    File vcf_in
    File? vcf_simu
    String program
    File reference
    File reference_idx
    File reference_dict
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
      vcf_norm_tbi = BiallelicNormalization.vcf_norm_tbi,
      vcf_simu = vcf_simu,
      reference = reference,
      reference_idx = reference_idx,
      reference_dict = reference_dict
  }

  output {
    File vcf_norm = BiallelicNormalization.vcf_norm
    File vcf_norm_tbi = BiallelicNormalization.vcf_norm_tbi
    File vcfEval = VariantEval.vcfEval
  }
}

# Split all markers by row in VCF
# Fix indels positions
task BiallelicNormalization {
  input {
    File vcf_file
    File reference
    File reference_idx
  }

  Int disk_size = ceil(size(vcf_file, "GB") + size(reference, "GB") + 2)

  command <<<
    bcftools norm ~{vcf_file} --rm-dup all -Ov --check-ref w -f ~{reference} > vcf_norm.vcf

    bgzip vcf_norm.vcf
    tabix -p vcf vcf_norm.vcf.gz
  >>>

  runtime {
    docker: "lifebitai/bcftools:1.10.2"
    #memory: "1 GB"
    #preemptible: 3
    #cpu: 1
    #disks: "local-disk " + disk_size + " HDD"
    job_name: "BiallelicNormalization"
    node:"--nodes=1"
    mem:"--mem=10G"
    tasks:"--ntasks=1"
    time:"01:00:00"
  }

  meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Uses [bcftools](https://samtools.github.io/bcftools/bcftools.html) to left-align and normalize indels in the VCF file."
  }

  output {
    File vcf_norm = "vcf_norm.vcf.gz"
    File vcf_norm_tbi = "vcf_norm.vcf.gz.tbi"
  }
}

task VariantEval {
  input {
    File vcf_norm
    File vcf_norm_tbi
    File? vcf_simu
    File reference
    File reference_idx
    File reference_dict
  }

  Int disk_size = ceil(size(vcf_norm, "GB") + size(reference, "GB") + size(vcf_simu, "GB") + 2)

  command <<<
    java -jar  /usr/gitc/GATK35.jar -T VariantEval -R ~{reference} -eval ~{vcf_norm} ~{"-D " + vcf_simu} -EV ValidationReport -EV CountVariants -o vcfEval.txt
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    #memory: "1 GB"
    #preemptible: 3
    #cpu: 1
    #disks: "local-disk " + disk_size + " HDD"
    job_name: "VariantEval"
    node:"--nodes=1"
    mem:"--mem=10G"
    tasks:"--ntasks=1"
    time:"01:00:00"
  }

  meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Uses [VariantEval](https://gatk.broadinstitute.org/hc/en-us/articles/360040507171-VariantEval-BETA-#:~:text=Overview,of%20s%20per%20sample%3B%20etc.) to generate report comparing variants estimated and simulated."
  }

  output {
    File vcfEval = "vcfEval.txt"
  }

}
