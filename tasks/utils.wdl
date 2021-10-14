version 1.0

import "../structs/struct_alignment.wdl"

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
    # memory: "2 GB"
    # cpu:1
    # preemptible: 3
    job_name: "ApplyRandomFilters"
    node:"--nodes=1"
    mem:"--mem=5GB"
    tasks:"--ntasks=1"
    time:"01:00:00"
  }

  output {
    File gatk_vcf_filt = "gatk_vcf_filt.vcf"
    File freebayes_vcf_filt = "freebayes_vcf_filt.vcf"
    File gatk_vcf_bam_counts_filt = "gatk_vcf_bam_counts_filt.vcf"
    File freebayes_vcf_bam_counts_filt = "freebayes_vcf_bam_counts_filt.vcf"
  }
}

task SplitMarkers {
  input{
    File vcf_file
  }

  command <<<
    bcftools view --max-alleles 2 --min-alleles 2 --output-type z --output-file biallelics.vcf.gz  ~{vcf_file}
    bcftools view --min-alleles 3 --types mnps --output-type z --output-file multiallelics.vcf.gz  ~{vcf_file}
  >>>

  runtime {
    docker:"lifebitai/bcftools:1.10.2"
    # memory: "2 GB"
    # cpu:1
    # preemptible: 3
    job_name: "SplitMarkers"
    node:"--nodes=1"
    mem:"--mem=5GB"
    tasks:"--ntasks=1"
    time:"01:00:00"
  }

  output {
    File biallelics = "biallelics.vcf.gz"
    File multiallelics = "multiallelics.vcf.gz"
  }
}

task JointMarkers{
  input{
    File biallelic_vcf
    File? multiallelic_vcf
  }

  command <<<

    filename=$(basename -- "~{biallelic_vcf}")
    extension="${filename##*.}"

    if [ "$extension" = "gz"]
    then 
      tabix -p vcf ~{biallelic_vcf}
      bcftools query -l ~{biallelic_vcf} | sort > samples.txt
      bcftools view -S samples.txt ~{biallelic_vcf} > biallelic_sort.vcf
    else 
      bgzip ~{biallelic_vcf}
      tabix -p vcf ~{biallelic_vcf}.gz
      bcftools query -l ~{biallelic_vcf}.gz | sort > samples.txt
      bcftools view -S samples.txt ~{biallelic_vcf}.gz > biallelic_sort.vcf
    fi
    
    tabix -p vcf ~{multiallelic_vcf}
    bcftools view -S samples.txt ~{multiallelic_vcf} > multiallelic_sort.vcf
    bcftools concat biallelic_sort.vcf multiallelic_sort.vcf --output merged.vcf
    bgzip merged.vcf

  >>>

  runtime {
    docker:"lifebitai/bcftools:1.10.2"
    # memory: "2 GB"
    # cpu:1
    # preemptible: 3
    job_name: "JointMarkers"
    node:"--nodes=1"
    mem:"--mem=5GB"
    tasks:"--ntasks=1"
    time:"01:00:00"
  }

  output {
    File merged_vcf = "merged.vcf.gz"
  }
}

# Only replace for biallelic markers
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

    bcftools view --min-alleles 3 ~{vcf} -Oz -o multiallelics.vcf.gz
    bcftools view -G --max-alleles 2 -v snps ~{vcf} -Oz -o sites.vcf.gz
    bcftools index --tbi -f sites.vcf.gz
    bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' sites.vcf.gz | bgzip -c > sites.tsv.gz
    bcftools mpileup -f ~{ref_fasta} -I -E -a 'FORMAT/DP,FORMAT/AD' -T sites.vcf.gz ~{sep=" " bams} -Ou > temp
    bcftools call temp -Aim -C alleles -T sites.tsv.gz  -o bam_vcf.vcf

    bcftools query -l multiallelics.vcf.gz | sort > samples.txt
    bcftools view -S samples.txt bam_vcf.vcf > biallelic_sort.vcf
    bcftools view -S samples.txt  multiallelics.vcf.gz > multiallelic_sort.vcf

    bgzip biallelic_sort.vcf
    tabix -p vcf biallelic_sort.vcf.gz

    bgzip multiallelic_sort.vcf
    tabix -p vcf multiallelic_sort.vcf.gz

    bcftools concat biallelic_sort.vcf.gz multiallelic_sort.vcf.gz -a -Oz --output ~{program}_bam_vcf.vcf.gz
    
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

task Compress {
    input{
      String name 
      Array[File] RDatas
      Array[File] maps_report
      Array[File] times
      Array[File] filters_report
      Array[File] errors_report
    }

    command <<<

      mkdir ~{name}

      mv ~{sep=" " RDatas} ~{sep=" " maps_report} \
        ~{sep=" " times} ~{sep=" " filters_report} \
        ~{sep=" " errors_report}  ~{name}

      tar -czvf ~{name}.tar.gz ~{name}

    >>>

  runtime{
    docker:"ubuntu:20.04"
    # memory: "2 GB"
    # cpu:1
    # preemptible: 3
    job_name: "Compress"
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
    input{
      String name 
      Array[File] RDatas
      Array[File] maps_report
      Array[File] times
    }
    
    command <<<

      mkdir ~{name}
      mv ~{sep=" " RDatas} ~{sep=" " maps_report} \
                ~{sep=" " times} ~{name}

      tar -czvf ~{name}.tar.gz ~{name}

    >>>

  runtime{
    docker:"ubuntu:20.04"
    # memory: "2 GB"
    # cpu:1
    # preemptible: 3
    job_name: "CompressGusmap"
    node:"--nodes=1"
    mem:"--mem=10GB"
    tasks:"--ntasks=1"
    time:"01:00:00"
  }

  output {
    File tar_gz_report = "~{name}.tar.gz"
  }

}


task GetMarkersPos {
  input{
    File true_vcf
    File filtered_gatk_vcf
    File filtered_gatk_vcf_bamcounts
    File filtered_freebayes_vcf
    File filtered_freebayes_vcf_bamcounts
  }

  command <<<

    bcftools query -f '%POS\n' ~{true_vcf} > true_vcf.tsv
    bcftools query -f '%POS\n' ~{filtered_gatk_vcf} > gatk_vcf_pos.tsv
    bcftools query -f '%POS\n' ~{filtered_gatk_vcf_bamcounts} > gatk_bam_pos.tsv
    bcftools query -f '%POS\n' ~{filtered_freebayes_vcf} > freebayes_vcf_pos.tsv
    bcftools query -f '%POS\n' ~{filtered_freebayes_vcf_bamcounts} > freebaye_bam_pos.tsv

    mkdir positions
    mv *tsv positions
    tar -czvf positions.tar.gz positions/
  >>>

  runtime { 
    docker:"lifebitai/bcftools:1.10.2"
    # memory: "2 GB"
    # cpu:1
    # preemptible: 3
    job_name: "GetMarkerPos"
    node:"--nodes=1"
    mem:"--mem=10GB"
    tasks:"--ntasks=1"
    time:"01:00:00"
  }

  output {
    File positions = "positions.tar.gz"
  }
}