version 1.0

import "../structs/alignment_struct.wdl"

task VcftoolsMerge {

  input {
    String prefix
    Array[File] vcfs
    Array[File] tbis
  }

  command <<<
    echo "~{sep=' ' tbis}"
    vcf-merge ~{sep=" "  vcfs} > ~{prefix}.variants.vcf
    bgzip ~{prefix}.variants.vcf
    tabix -p vcf ~{prefix}.variants.vcf.gz

  >>>
  runtime {
    docker: "cristaniguti/split_markers:0.0.1"
    memory:"4 GB"
    cpu:1
    preemptible: 3
  }

  output {
    File vcf = "~{prefix}.variants.vcf.gz"
    File tbi = "~{prefix}.variants.vcf.gz.tbi"
  }
}

task VcftoolsApplyFilters {

  input {
    File vcf_in
    Float max_missing
    Int min_alleles
    Int max_alleles
    Float? maf
    String program
    Int? min_meanDP
    String? chromosome
  }

  command <<<
    vcftools --gzvcf ~{vcf_in} --max-missing ~{max_missing} --min-alleles ~{min_alleles} --max-alleles ~{max_alleles} ~{"--maf " +  maf} ~{"--min-meanDP " +  min_meanDP} ~{"--chr " +  chromosome} --recode --out ~{program}

    bgzip ~{program}.recode.vcf
    tabix -p vcf ~{program}.recode.vcf.gz

  >>>
  runtime {
    docker: "cristaniguti/split_markers:0.0.1"
    memory:"5 GB"
    cpu:1
    preemptible: 3
  }

  output {
    File vcf = "~{program}.recode.vcf.gz"
    File tbi = "~{program}.recode.vcf.gz.tbi"
  }
}

task ApplyRandomFilters {
  input{
    File gatk_vcf
    File freebayes_vcf
    File gatk_vcf_bam_counts
    File freebayes_vcf_bam_counts
    String? filters
    String? chromosome
  }

  command <<<
    vcftools --gzvcf ~{gatk_vcf} ~{filters} ~{"--chr " + chromosome} --recode --stdout > gatk_vcf_filt.vcf

    vcftools --gzvcf ~{freebayes_vcf} ~{filters} ~{"--chr " + chromosome} --recode --stdout > freebayes_vcf_filt.vcf

    vcftools --gzvcf ~{gatk_vcf_bam_counts} ~{filters} ~{"--chr " + chromosome} --recode --stdout > gatk_vcf_bam_counts_filt.vcf

    vcftools --gzvcf ~{freebayes_vcf_bam_counts} ~{filters} ~{"--chr " + chromosome} --recode --stdout > freebayes_vcf_bam_counts_filt.vcf
  >>>

  runtime {
    docker:"cristaniguti/split_markers:0.0.1"
    memory: "2 GB"
    cpu:1
    preemptible: 3
  }

  output {
    File gatk_vcf_filt = "gatk_vcf_filt.vcf"
    File freebayes_vcf_filt = "freebayes_vcf_filt.vcf"
    File gatk_vcf_bam_counts_filt = "gatk_vcf_bam_counts_filt.vcf"
    File freebayes_vcf_bam_counts_filt = "freebayes_vcf_bam_counts_filt.vcf"
  }
}

task ReplaceAD {
  input {
    File ref_fasta
    File ref_index
    Array[File] bams
    Array[File] bais
    File vcf
    File tbi
    String program
  }

  command <<<

    bcftools view -G -v snps ~{vcf} -Oz -o sites.vcf.gz
    bcftools index --tbi -f sites.vcf.gz
    bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' sites.vcf.gz | bgzip -c > sites.tsv.gz
    bcftools mpileup -f ~{ref_fasta} -I -E -a 'FORMAT/DP,FORMAT/AD' -T sites.vcf.gz ~{sep=" " bams} -Ou > temp
    bcftools call temp -Aim -C alleles -T sites.tsv.gz  -o ~{program}_bam_vcf.vcf

    bgzip ~{program}_bam_vcf.vcf
  >>>

  runtime {
    docker:"lifebitai/bcftools:1.10.2"
    # memory: "2 GB"
    # cpu:1
    # preemptible: 3
    job_name: "ReplaceAD"
    node:"--nodes=1"
    mem:"--mem=50GB"
    tasks:"--ntasks=1"
    time:"24:00:00"
  }

  output {
    File bam_vcf =  "~{program}_bam_vcf.vcf.gz"
  }
}


# Deprecated
task MergeVCFLines {
  input {
    File vcf_file
    File reference
    File reference_idx
  }

  command <<<

    bcftools norm -m+any ~{vcf_file} -f ~{reference} -o multiallelics.vcf.gz

  >>>

  runtime {
    docker:"lifebitai/bcftools:1.10.2"
    # memory: "2 GB"
    # cpu:1
    # preemptible: 3
    job_name: "ReplaceAD"
    node:"--nodes=1"
    mem:"--mem=50GB"
    tasks:"--ntasks=1"
    time:"24:00:00"
  }

  output {
    File multi_vcf = "multiallelics.vcf.gz"
  }
}


# Phase VCF file according with WhatsHap
# Require BAM files
task RunWhatsHap {
  input{
    File merged_bam
    File reference
    File vcf_file
  }

  command <<<
    whatshap phase -o phased.vcf --reference=~{reference} ~{vcf_file} ~{merged_bam} --tag=PS
  >>>

  runtime{
    docker:"cristaniguti/whatshap:0.0.1"
    # memory: "2 GB"
    # cpu:1
    # preemptible: 3
    job_name: "WhatsHap"
    node:"--nodes=1"
    mem:"--mem=50GB"
    tasks:"--ntasks=1"
    time:"24:00:00"
  }

  output{
    phased_vcf = "phased.vcf"
  }
}


task Compress {
    input:
      String name 
      Array[File] RDatas
      Array[File] maps_report
      Array[File] times
      Array[File] filters_report
      Array[File] errors_report

    command <<<

      tar -czvf ~{name}.tar.gz ~{sep=" " RDatas} ~{sep=" " maps_report} \
                ~{sep=" " times} ~{sep=" " filters_report} ~{sep=" " errors_report}

    >>>

  runtime{
    docker:"cristaniguti/ubuntu:20.04"
    # memory: "2 GB"
    # cpu:1
    # preemptible: 3
    job_name: "WhatsHap"
    node:"--nodes=1"
    mem:"--mem=10GB"
    tasks:"--ntasks=1"
    time:"01:00:00"
  }

  output {
    File tar_gz_report = "~{name}.tar.gz"
  }

}


task CompressGusmap {
    input:
      String name 
      Array[File] RDatas
      Array[File] maps_report
      Array[File] times

    command <<<

      tar -czvf ~{name}.tar.gz ~{sep=" " RDatas} ~{sep=" " maps_report} \
                ~{sep=" " times} 

    >>>

  runtime{
    docker:"cristaniguti/ubuntu:20.04"
    # memory: "2 GB"
    # cpu:1
    # preemptible: 3
    job_name: "WhatsHap"
    node:"--nodes=1"
    mem:"--mem=10GB"
    tasks:"--ntasks=1"
    time:"01:00:00"
  }

  output {
    File tar_gz_report = "~{name}.tar.gz"
  }

}