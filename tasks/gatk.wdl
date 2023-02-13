version 1.0

## Process all samples because it RAD experiments
## usually do not have large ammount of reads.
## NotE: if BAMS have same name it will be overrided in this task.
task HaplotypeCaller {
  input {
    File reference_fasta
    File reference_dict
    File reference_fai
    Array[File] bams
    Array[File] bams_index
    Int ploidy
    Int chunk_size
  }

  Int disk_size = ceil((size(bams, "GiB") + 30) + size(reference_fasta, "GiB")) + 20
  Int memory_max = ceil(5000 * chunk_size)
  Int memory_min = memory_max / 2
  Int memory_size = memory_max + 5000
  Int max_cores = ceil(chunk_size * 4 + 2)

  command <<<
    set -euo pipefail

    for bam in ~{sep=" " bams}; do ln -s $bam .; done
    for bai in ~{sep=" " bams_index}; do ln -s $bai .; done

    mkdir vcfs
    ## gvcf for each sample
    for bam in *.bam; do
      out_name=$(basename -s ".bam" "$bam")
      /usr/gitc/gatk4/./gatk --java-options "-Xms~{memory_min}m -Xmx~{memory_max}m" HaplotypeCaller \
        -ERC GVCF \
        -R ~{reference_fasta} \
        -ploidy ~{ploidy} \
        -I "$bam" \
        -O "vcfs/${out_name}.g.vcf.gz" \
        --max-alternate-alleles 1 \
        --max-reads-per-alignment-start 0 &
    done

    wait  # TODO: Why this line? Because of the &
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    cpu: max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "HaplotypeCaller"
    mem:"~{memory_size}M"
    time:"24:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Runs [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)."
  }

  output {
    Array[File] vcfs = glob("vcfs/*.vcf.gz")
    Array[File] vcfs_index = glob("vcfs/*.vcf.gz.tbi")
  }
}

task ImportGVCFs  {
  input {
    Array[File] vcfs
    Array[File] vcfs_index
    File reference_fasta
    File reference_fai
    File reference_dict
    String interval
  }

  Int disk_size = ceil(size(vcfs, "GiB") * 1.5 + size(reference_fasta, "GiB") * 1.5)
  Int memory_max = 2300
  Int memory_min = 2000
  Int memory_size = 26000

  command <<<
    set -euo pipefail
    grep ">" ~{reference_fasta} | sed 's/^.//' > interval.list
    mkdir gvcfs
    for i in ~{sep=" " vcfs}; do ln -s $i gvcfs/; done
    for i in ~{sep=" " vcfs_index}; do ln -s $i gvcfs/; done

    /usr/gitc/gatk4/./gatk --java-options "-Xms~{memory_min}m -Xmx~{memory_max}m" GenomicsDBImport \
      --batch-size 50 \
      --reader-threads 5 \
      --genomicsdb-workspace-path cohort_db \
      -L ~{interval} \
      -V $(find gvcfs/*.g.vcf.gz -type l | paste -d',' -s | sed 's/,/ -V /g') \
      --consolidate 

    tar -cf cohort_db.tar cohort_db

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    cpu: 4
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "ImportGVCFs"
    mem:"~{memory_size}M"
    time:"24:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Runs [GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport)."
  }

  output {
    File output_workspace = "cohort_db.tar"
  }
}

task GenotypeGVCFs   {
  input {
    File workspace_tar
    File reference_fasta
    File reference_fai
    File reference_dict
    String interval
  }

  Int disk_size = ceil(size(reference_fasta, "GiB") * 1.5 + size(workspace_tar, "GiB") * 1.5)
  Int memory_max = 2300
  Int memory_min = 2000
  Int memory_size = 26000

  command <<<
    set -euo pipefail

    tar -xf ~{workspace_tar}

    /usr/gitc/gatk4/./gatk --java-options "-Xms~{memory_min}m -Xmx~{memory_max}m" GenotypeGVCFs \
      -R ~{reference_fasta} \
      -V gendb://cohort_db \
      -L ~{interval} \
      -G StandardAnnotation \
      --max-alternate-alleles 1 \
      -O gatk.vcf.gz 

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    cpu: 2
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "GenotypeGVCFs"
    mem:"~{memory_size}M"
    time:"24:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Runs [GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs)."
  }

  output {
    File vcf = "gatk.vcf.gz"
    File vcf_tbi = "gatk.vcf.gz.tbi"
  }
}

task MergeVCFs {

  input {
    Array[File] input_vcfs
    Array[File] input_vcf_indices
  }

  Int disk_size = ceil(size(input_vcfs, "GiB") * 2.5) + 10
  Int memory_max = 2800
  Int memory_min = 2600
  Int memory_size = 3000

  command <<<

    /usr/gitc/gatk4/./gatk --java-options "-Xms~{memory_min}m -Xmx~{memory_max}m" \
      MergeVcfs \
        -I ~{sep=' -I' input_vcfs} \
        -O gatk_joint.vcf.gz
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    cpu: 1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "MergeVCFs"
    mem:"~{memory_size}M"
    time:"10:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Runs [MergeVcfs](https://gatk.broadinstitute.org/hc/en-us/articles/360037226612-MergeVcfs-Picard-)."
  }

  output {
    File output_vcf = "gatk_joint.vcf.gz"
    File output_vcf_index = "gatk_joint.vcf.gz.tbi"
  }
}

task VariantsToTable {
    input {
        File vcf_file
        File vcf_tbi
        File reference
        File reference_dict
        File reference_idx
    }

    Int disk_size = ceil(size(reference, "GB") + size(vcf_file, "GB") + 2)
    Int memory_size = 2500

    command <<<

        /usr/gitc/gatk4/./gatk VariantsToTable \
            -V ~{vcf_file} \
            -F CHROM -F POS \
            -F QD \
            -F FS \
            -F SOR \
            -F MQ \
            -F MQRankSum \
            -F ReadPosRankSum \
            -O Total.table

    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
        cpu: 1
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "VariantsToTable"
        mem:"~{memory_size}M"
        time:"01:40:00"
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Generated tables with estimated markers quality parameters. See more information in [VariatsToTable](https://gatk.broadinstitute.org/hc/en-us/articles/360036896892-VariantsToTable) tool"
    }

    output {
        File Total = "Total.table"
    }

}

task VariantFiltration {
    input {
        File vcf_file
        File vcf_tbi
        File reference
        File reference_idx
        File reference_dict
    }

    Int disk_size = ceil(size(vcf_file, "GB") + size(reference, "GB") + 1)
    Int memory_size = 3000

    command <<<
        /usr/gitc/gatk4/./gatk VariantFiltration \
            -V ~{vcf_file} \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O gatk_filters.vcf.gz

        /usr/gitc/gatk4/./gatk SelectVariants \
            -R ~{reference} \
            -V gatk_filters.vcf.gz \
            --exclude-filtered \
            -O gatk_filtered.vcf.gz

    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
        cpu: 1
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "VariantFiltration"
        mem:"~{memory_size}M"
        time:"01:00:00"
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Filters empirical VCF according to GATK Hard Filtering. See more information in [VariatsToTable](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants) tool"
    }

    output {
        File filters_vcf = "gatk_filters.vcf.gz"
        File filtered_vcf = "gatk_filtered.vcf.gz"
    }
}

task VariantsToTableForHardFilteringSimulated {
    input {
        File vcf_file
        File vcf_tbi
        File? simu_vcf
        File reference
        File reference_dict
        File reference_idx
    }

     Int disk_size = ceil(size(reference, "GB") + size(vcf_file, "GB") + size(simu_vcf, "GB") + 2)
     Int memory_size = 3000

    command <<<
        /usr/gitc/./bgzip -c  ~{simu_vcf} > ~{simu_vcf}.gz
        /usr/gitc/./tabix -p vcf ~{simu_vcf}.gz

        /usr/gitc/gatk4/./gatk SelectVariants \
            -R ~{reference} \
            -V ~{vcf_file} \
            --discordance ~{simu_vcf}.gz \
            -O FalsePositives.vcf.gz

        /usr/gitc/gatk4/./gatk SelectVariants \
            -R ~{reference} \
            -V ~{vcf_file} \
            --concordance ~{simu_vcf}.gz \
            -O TruePositives.vcf.gz

        /usr/gitc/gatk4/./gatk VariantsToTable \
            -V ~{vcf_file} \
            -F CHROM -F POS \
            -F QD \
            -F FS \
            -F SOR \
            -F MQ \
            -F MQRankSum \
            -F ReadPosRankSum \
            -O Total.table

        /usr/gitc/gatk4/./gatk VariantsToTable \
            -V FalsePositives.vcf.gz \
            -F CHROM -F POS \
            -F QD \
            -F FS \
            -F SOR \
            -F MQ \
            -F MQRankSum \
            -F ReadPosRankSum \
            -O FalsePositives.table

        /usr/gitc/gatk4/./gatk VariantsToTable \
            -V TruePositives.vcf.gz \
            -F CHROM -F POS \
            -F QD \
            -F FS \
            -F SOR \
            -F MQ \
            -F MQRankSum \
            -F ReadPosRankSum \
            -O TruePositives.table
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
        cpu: 1
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "VariantsToTable"
        mem:"~{memory_size}M"
        time:"01:40:00"
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Generated tables with simulated and estimated markers quality parameters. See more information in [VariatsToTable](https://gatk.broadinstitute.org/hc/en-us/articles/360036896892-VariantsToTable) tool"
    }

    output {
        File FalsePositives = "FalsePositives.table"
        File TruePositives = "TruePositives.table"
        File Total = "Total.table"
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
    Int ploidy
  }

  Int disk_size = ceil(size(vcf_norm, "GiB") + size(reference, "GiB") + size(vcf_simu, "GiB") + 2)
  Int memory_size = 3000

  command <<<
    java -jar  /usr/gitc/GATK35.jar -T VariantEval -R ~{reference} -eval ~{vcf_norm} ~{"-D " + vcf_simu} -EV ValidationReport -EV CountVariants -ploidy ~{ploidy} -o vcfEval.txt
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    cpu: 1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "VariantEval"
    mem:"~{memory_size}M"
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
