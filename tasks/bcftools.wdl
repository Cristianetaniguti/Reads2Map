version 1.0


# Split all markers by row in VCF
# Fix indels positions
task BiallelicNormalization {
  input {
    File vcf_file
    File reference
    File reference_idx
  }

  Int disk_size = ceil(size(vcf_file, "GiB") + size(reference, "GiB") + 2)
  Int memory_size = 4000 + ceil(size(vcf_file, "MiB") + size(reference, "MiB") + 2)

  command <<<

    bcftools norm ~{vcf_file} --rm-dup all -Ov --check-ref w -f ~{reference} > vcf_norm.vcf

    bgzip vcf_norm.vcf
    tabix -p vcf vcf_norm.vcf.gz
  >>>

  runtime {
    docker: "lifebitai/bcftools:1.10.2"
    singularity:"docker://lifebitai/bcftools:1.10.2"
    cpu: 1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "BiallelicNormalization"
    mem:"~{memory_size}M"
    time: 1
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


task FixTasselVCF {
  input {
    File vcf_file
    File reference
    File reference_idx
  }

  Int disk_size = ceil(size(vcf_file, "GiB") + size(reference, "GiB") + 2)
  Int memory_size = 4000 + ceil(size(vcf_file, "MiB") + size(reference, "MiB") + 2)

  command <<<

    sed 's/PL,Number=./PL,Number=G/g' ~{vcf_file} > tassel_fix.vcf
    sed -i 's/AD,Number=./AD,Number=R/g' tassel_fix.vcf
    sed -i 's/AF,Number=./AF,Number=A/g' tassel_fix.vcf 
    sed -i '/INFO=<ID=AF/a ##INFO=<ID=QualityScore,Number=1,Type=Float,Description="Tassel specific score">' tassel_fix.vcf 

    grep ">" ~{reference} > chrs
    sed -i 's/>//' chrs
    cp chrs chrs_tassel
    sed -i 's/Chr//' chrs_tassel
    sed -i 's/scaffold/SCAFFOLD/' chrs_tassel
    paste -d'\t' chrs_tassel chrs > fix_chrom

    bgzip tassel_fix.vcf
    tabix -p vcf tassel_fix.vcf.gz
    bcftools annotate --rename-chrs fix_chrom tassel_fix.vcf.gz > tassel_fix_chr.vcf.gz

  >>>

  runtime {
    docker: "lifebitai/bcftools:1.10.2"
    singularity:"docker://lifebitai/bcftools:1.10.2"
    cpu: 1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "BiallelicNormalization"
    mem:"~{memory_size}M"
    time: 1
  }

  meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Uses [bcftools](https://samtools.github.io/bcftools/bcftools.html) to left-align and normalize indels in the VCF file."
  }

  output {
    File vcf_fixed = "tassel_fix_chr.vcf.gz"
  }
}