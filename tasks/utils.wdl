version 1.0

task mergeVCFs {
    input {
        Array[File] haplo_vcf
    }

    Int disk_size = ceil(size(haplo_vcf, "GiB") * 2)
    Int memory_size =  ceil(3000 + size(haplo_vcf, "MiB") * 2)

    command <<<

        vcfs=(~{sep = " " haplo_vcf})

        index=1
        for file in ${!vcfs[*]}; do 
            filename=$(basename -- "${vcfs[$file]}")
            name="${filename%.*}_$index"
            echo $name
            tabix -p vcf ${vcfs[$file]} 
            bcftools sort ${vcfs[$file]} --output-file $name.sorted.vcf
            bgzip $name.sorted.vcf
            tabix -p vcf $name.sorted.vcf.gz
            index=$(expr $index + 1)
        done

        bcftools concat *sorted.vcf.gz --output merged.vcf
        bcftools sort merged.vcf --output-file merged.sorted.vcf
        bgzip merged.sorted.vcf

    >>>

    runtime {
        docker:"lifebitai/bcftools:1.10.2"
        singularity: "docker://lifebitai/bcftools:1.10.2"
        cpu: 1
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "mergeVCFs"
        mem:"~{memory_size}M"
        time: 1
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Uses [bcftools](https://samtools.github.io/bcftools/bcftools.html) to sort and merge VCF files."
    }

    output {
        File merged_vcf = "merged.sorted.vcf.gz"
    }
}

# It will always produce P1, P2, F1 and then F2_00X, where
# X will increase from 1 to samples
task GenerateSampleNames {  # TODO: probably a name like 'ReadSamplesNamesInVcf' is better

  input {
    File simulated_vcf
  }

  Int disk_size = ceil(size(simulated_vcf, "GiB") * 2) 
  Int memory_size = ceil(1000 + size(simulated_vcf, "MiB") * 2)  

  command <<<
    export PATH=$PATH:/opt/conda/bin

    python <<CODE
    from pysam import VariantFile

    bcf_in = VariantFile("~{simulated_vcf}")

    for i in bcf_in.header.samples:
        print(i)
    CODE

  >>>

  runtime {
    docker: "cristaniguti/miniconda-alpine:0.0.1"
    singularity: "docker://cristaniguti/miniconda-alpine:0.0.1"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "GenerateSampleNames"
    mem:"~{memory_size}M"
    time: 5
  }

  meta {
    author: "Lucas Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Creates the sample names."
  }

  output {
    Array[String] names = read_lines(stdout())
  }
}

task ApplyRandomFilters {
  input{
    File? gatk_vcf
    File? freebayes_vcf
    File? gatk_vcf_bam_counts
    File? freebayes_vcf_bam_counts
    String? filters
    String? chromosome
  }

  Int disk_size = ceil(size(gatk_vcf, "GiB") * 2 + size(freebayes_vcf, "GiB") * 2 + size(gatk_vcf_bam_counts, "GiB") * 2 + size(freebayes_vcf_bam_counts, "GiB") * 2)
  Int memory_size = ceil(3000 + size(gatk_vcf, "MiB") * 2 + size(freebayes_vcf, "MiB") * 2 + size(gatk_vcf_bam_counts, "MiB") * 2 + size(freebayes_vcf_bam_counts, "MiB") * 2)

  command <<<
    # Required update to deal with polyploids
    zcat  ~{gatk_vcf} | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' > out1.vcf
    zcat ~{gatk_vcf_bam_counts} | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/'  > out2.vcf

    vcftools --gzvcf out1.vcf  ~{filters} ~{"--chr " + chromosome} --recode --stdout > gatk_vcf_filt.vcf
    vcftools --gzvcf out2.vcf ~{filters} ~{"--chr " + chromosome} --recode --stdout > gatk_vcf_bam_counts_filt.vcf

    vcftools --gzvcf ~{freebayes_vcf} ~{filters} ~{"--chr " + chromosome} --recode --stdout > freebayes_vcf_filt.vcf
    vcftools --gzvcf ~{freebayes_vcf_bam_counts} ~{filters} ~{"--chr " + chromosome} --recode --stdout > freebayes_vcf_bam_counts_filt.vcf
  >>>

  runtime {
    docker:"cristaniguti/split_markers:0.0.1"
    singularity: "docker://cristaniguti/split_markers:0.0.1"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "ApplyRandomFilters"
    mem:"~{memory_size}M"
    time: 1
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Uses [vcftools](http://vcftools.sourceforge.net/) to filter VCF file by user-defined criterias."
  }

  output {
    File gatk_vcf_filt = "gatk_vcf_filt.vcf"
    File freebayes_vcf_filt = "freebayes_vcf_filt.vcf"
    File gatk_vcf_bam_counts_filt = "gatk_vcf_bam_counts_filt.vcf"
    File freebayes_vcf_bam_counts_filt = "freebayes_vcf_bam_counts_filt.vcf"
  }
}

task ApplyRandomFiltersArray {
  input{
    Array[File] vcfs
    Array[String] vcfs_SNPCall_software
    Array[String] vcfs_Counts_source
    Array[String] vcfs_GenoCall_software
    String? filters = "-r " + chromosome 
    String? chromosome
  }

  Int disk_size = ceil(size(vcfs, "GiB") * 2)
  Int memory_size = ceil(3000 + size(vcfs, "MiB") * 2)

  command <<<

      vcfs=(~{sep = " " vcfs})
      vcfs_snp_software=(~{sep=" " vcfs_SNPCall_software})
      vcfs_counts_source=(~{sep=" " vcfs_Counts_source})
      vcfs_geno_software=(~{sep=" " vcfs_GenoCall_software})

      for index in ${!vcfs[*]}; do
          if [[ ${vcfs[$index]} != *.gz ]]; then
            cp ${vcfs[$index]} temp.vcf
            bgzip temp.vcf
          else 
            cp ${vcfs[$index]} temp.vcf.gz
          fi
          
          tabix -p vcf temp.vcf.gz
          bcftools view temp.vcf.gz ~{filters} \
          -o vcf_filt_${vcfs_snp_software[$index]}_${vcfs_counts_source[$index]}_${vcfs_geno_software[$index]}.vcf
          bgzip vcf_filt_${vcfs_snp_software[$index]}_${vcfs_counts_source[$index]}_${vcfs_geno_software[$index]}.vcf
          rm temp.vcf.gz temp.vcf.gz.tbi
          echo vcf_filt_${vcfs_snp_software[$index]}_${vcfs_counts_source[$index]}_${vcfs_geno_software[$index]}.vcf.gz >> outputs.txt
      done

  >>>

  runtime {
    docker:"lifebitai/bcftools:1.10.2"
    singularity: "docker://lifebitai/bcftools:1.10.2"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "ApplyRandomFilters"
    mem:"~{memory_size}M"
    time: 1
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Uses [vcftools](http://vcftools.sourceforge.net/) to filter VCF file by user-defined criterias."
  }

  output {
    Array[File] vcfs_filt = read_lines("outputs.txt")
  }
}


task SplitMarkers {
  input{
    File vcf_file
  }

  Int disk_size = ceil(size(vcf_file, "GiB") * 2)
  Int memory_size = ceil(3000 + size(vcf_file, "MiB") * 2)

  command <<<
    bcftools view --max-alleles 2 --min-alleles 2 --output-type z --output-file biallelics.vcf.gz  ~{vcf_file}
    bcftools view --min-alleles 3 --types mnps --output-type z --output-file multiallelics.vcf.gz  ~{vcf_file}
  >>>

  runtime {
    docker:"lifebitai/bcftools:1.10.2"
    singularity: "docker://lifebitai/bcftools:1.10.2"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "SplitMarkers"
    mem:"~{memory_size}M"
    time: 1
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Uses [bcftools](https://samtools.github.io/bcftools/bcftools.html) to split the VCF in biallelic and multiallelic markers."
  }

  output {
    File biallelics = "biallelics.vcf.gz"
    File multiallelics = "multiallelics.vcf.gz"
  }
}

task JointMarkers {
  input{
    String SNPCall_program
    String CountsFrom
    String GenotypeCall_program
    File biallelic_vcf
    File? multiallelic_vcf
  }

   Int disk_size = ceil(size(biallelic_vcf, "GiB") * 2 + size(multiallelic_vcf, "GiB") * 2)
  Int memory_size = ceil(3000 + size(biallelic_vcf, "MiB") * 2 + size(multiallelic_vcf, "MiB") * 2)

  command <<<

    filename=$(basename -- "~{biallelic_vcf}")
    extension="${filename##*.}"

    if [ "$extension" = "gz" ]
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
    bcftools concat biallelic_sort.vcf multiallelic_sort.vcf --output ~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.vcf
    bgzip ~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.vcf

  >>>

  runtime {
    docker:"lifebitai/bcftools:1.10.2"
    singularity: "docker://lifebitai/bcftools:1.10.2"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "JointMarkers"
    mem:"~{memory_size}M"
    time: 1
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Uses [bcftools](https://samtools.github.io/bcftools/bcftools.html) to joint the VCF files with biallelic and multiallelic markers."
  }

  output {
    File merged_vcf = "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.vcf.gz"
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
    String counts_source
  }

  Int disk_size = ceil(size(ref_fasta, "GiB") + size(bams, "GiB") * 1.5 + size(vcf, "GiB") * 1.5)
  Int memory_size = ceil(size(vcf, "MiB") * 3 + 5000)

  command <<<

    bcftools view --min-alleles 3 ~{vcf} -Oz -o multiallelics.vcf.gz
    bcftools view -G --max-alleles 2 -v snps ~{vcf} -Oz -o sites.vcf.gz
    bcftools index --tbi -f sites.vcf.gz
    bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' sites.vcf.gz | bgzip -c > sites.tsv.gz
    bcftools mpileup -f ~{ref_fasta} -d 500000 -I -E -a 'FORMAT/DP,FORMAT/AD' -T sites.vcf.gz ~{sep=" " bams} -Ou > temp
    bcftools call temp -Aim -C alleles -T sites.tsv.gz  -o bam_vcf.vcf

    bcftools query -l multiallelics.vcf.gz | sort > samples.txt
    bcftools view -S samples.txt bam_vcf.vcf > biallelic_sort.vcf
    bcftools view -S samples.txt  multiallelics.vcf.gz > multiallelic_sort.vcf

    bgzip biallelic_sort.vcf
    tabix -p vcf biallelic_sort.vcf.gz

    bgzip multiallelic_sort.vcf
    tabix -p vcf multiallelic_sort.vcf.gz

    bcftools concat biallelic_sort.vcf.gz multiallelic_sort.vcf.gz -a -Oz --output ~{program}_bam_vcf.vcf.gz
    tabix -p vcf ~{program}_bam_vcf.vcf.gz

  >>>

  runtime {
    docker:"lifebitai/bcftools:1.10.2"
    singularity: "docker://lifebitai/bcftools:1.10.2"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "ReplaceAD"
    mem:"~{memory_size}M"
    time: 24
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Uses [bcftools](https://samtools.github.io/bcftools/bcftools.html) to replace Allele Depth VCF field with read depth information from BAM alignment files."
  }

  output {
    File bam_vcf =  "~{program}_bam_vcf.vcf.gz"
    File bam_vcf_tbi = "~{program}_bam_vcf.vcf.gz.tbi"
    String software = "~{program}"
    String source = "~{counts_source}"
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

    Int disk_size = ceil(size(RDatas, "GiB") + size(maps_report, "GiB") + size(times, "GiB") + size(filters_report, "GiB") + size(errors_report, "GiB"))
    Int memory_size = ceil(2000 + size(RDatas, "MiB") + size(maps_report, "MiB") + size(times, "MiB") + size(filters_report, "MiB") + size(errors_report, "MiB"))

    command <<<

      mkdir ~{name}

      cp ~{sep=" " RDatas} ~{sep=" " maps_report} \
        ~{sep=" " times} ~{sep=" " filters_report} \
        ~{sep=" " errors_report}  ~{name}

      tar -czvf ~{name}.tar.gz ~{name}

    >>>

  runtime {
    docker:"ubuntu:20.04"
    singularity: "docker://ubuntu:20.04"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "Compress"
    mem:"~{memory_size}M"
    time: 1
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Move resulted reports to a single directory and compress it."
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
    Int depth
    Int seed
  }

  Int disk_size = ceil(size(true_vcf, "GiB") * 1.5 + size(filtered_gatk_vcf, "GiB") * 1.5 + size(filtered_gatk_vcf_bamcounts, "GiB") + size(filtered_freebayes_vcf, "GiB") + size(filtered_freebayes_vcf_bamcounts, "GiB"))
  Int memory_size = ceil(3000 + size(true_vcf, "MiB") * 1.5 + size(filtered_gatk_vcf, "MiB") * 1.5 + size(filtered_gatk_vcf_bamcounts, "MiB") + size(filtered_freebayes_vcf, "MiB") + size(filtered_freebayes_vcf_bamcounts, "MiB"))

  command <<<

    bcftools query -f '%POS\n' ~{true_vcf} > ~{depth}_~{seed}_true_vcf.tsv
    bcftools query -f '%POS\n' ~{filtered_gatk_vcf} > ~{depth}_~{seed}_gatk_vcf_pos.tsv
    bcftools query -f '%POS\n' ~{filtered_gatk_vcf_bamcounts} > ~{depth}_~{seed}_gatk_bam_pos.tsv
    bcftools query -f '%POS\n' ~{filtered_freebayes_vcf} > ~{depth}_~{seed}_freebayes_vcf_pos.tsv
    bcftools query -f '%POS\n' ~{filtered_freebayes_vcf_bamcounts} > ~{depth}_~{seed}_freebaye_bam_pos.tsv

    mkdir ~{depth}_~{seed}_positions
    mv *tsv ~{depth}_~{seed}_positions
    tar -czvf ~{depth}_~{seed}_positions.tar.gz ~{depth}_~{seed}_positions/
  >>>

  runtime {
    docker:"lifebitai/bcftools:1.10.2"
    singularity: "docker://lifebitai/bcftools:1.10.2"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "GetMarkerPos"
    mem:"~{memory_size}M"
    time: 1
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Uses [bcftools](https://samtools.github.io/bcftools/bcftools.html) to extract markers position information."
  }

  output {
    File positions = "~{depth}_~{seed}_positions.tar.gz"
  }
}

task MergeBams{
    input {
        Array[File] bam_files
    }

    Int disk_size = ceil(size(bam_files, "GiB") * 2)
    Int memory_size = ceil(3000 + size(bam_files, "MiB") * 2)

    command <<<
        samtools merge merged.bam ~{sep=" " bam_files}
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
        singularity: "docker://us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
        cpu:1
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "MergeBams"
        mem:"~{memory_size}M"
        time: 10
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Uses [samtools](https://samtools.github.io/) to merge BAM alignment files."
    }

    output {
        File merged_bam = "merged.bam"
    }
}


task TarFiles {
  input {
    Array[File] sequences
  }

  Int disk_size = ceil(size(sequences, "GiB") * 2)
  Int memory_size = ceil(4000 + size(sequences, "MiB") * 2)
  
  command <<<
    mkdir results
    mv ~{sep=" " sequences} results
    tar -czvf results.tar.gz results
  >>>

  runtime {
    docker:"kfdrc/cutadapt"
    singularity: "docker://kfdrc/cutadapt"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "TarFiles"
    mem:"~{memory_size}M"
    time: 10
  }

  output {
    File results = "results.tar.gz"
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
    Int memory_size = ceil(5000 + size(vcf_file, "MiB") * 2)

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
        singularity: "docker://us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
        cpu: 1
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "VariantFiltration"
        mem:"~{memory_size}M"
        time: 1
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Filters simulated VCF according to GATK Hard Filtering. See more information in [VariatsToTable](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants) tool"
    }

    output {
        File filters_vcf = "gatk_filters.vcf.gz"
        File filtered_vcf = "gatk_filtered.vcf.gz"
    }
}

task BamToBed {
    input {
        File? merged_bams
    }

    Int disk_size = ceil(size(merged_bams, "GiB") * 1.5)
    Int memory_size = ceil(3000 + size(merged_bams, "MiB") * 1.5)

    command <<<
        bamToBed -i ~{merged_bams} > file.bed
        bedtools merge -i file.bed > merged.bed
    >>>

    runtime {
        docker: "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
        singularity: "docker://biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
        cpu: 1
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "BamToBed"
        mem:"~{memory_size}M"
        time: 5
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Uses [bamToBed](https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html) and [bedtools](https://bedtools.readthedocs.io/en/latest/) to create BED file and merge overlapping intervals."
    }

    output {
        File merged_bed = "merged.bed"
    }
}