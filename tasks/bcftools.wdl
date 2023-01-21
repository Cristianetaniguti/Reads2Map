version 1.0


# Split all markers by row in VCF
# Fix indels positions
task BiallelicNormalization {
  input {
    File vcf_file
    File reference
    File reference_idx
    Int ploidy
    String software
  }

  Int disk_size = ceil(size(vcf_file, "GiB") + size(reference, "GiB") + 2)
  Int memory_size = 7000

  command <<<

    #if [ ~{ploidy} -gt 2 ] && [ "~{software}" == "freebayes" ] # GATK returns an error when trying to split by row saying that PL has wrong number of fields
    #then
    #  bcftools norm ~{vcf_file} -m - -Ov --check-ref w -f ~{reference} > vcf_norm.vcf
    #else
      bcftools norm ~{vcf_file} --rm-dup all -Ov --check-ref w -f ~{reference} > vcf_norm.vcf
    #fi

    bgzip vcf_norm.vcf
    tabix -p vcf vcf_norm.vcf.gz
  >>>

  runtime {
    docker: "lifebitai/bcftools:1.10.2"
    cpu: 1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "BiallelicNormalization"
    mem:"~{memory_size}M"
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
