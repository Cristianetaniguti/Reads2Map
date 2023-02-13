version 1.0

import "../structs/dna_seq_structs.wdl"

# This task considers that it is possible to receive more than one fastq file per sample
# It keeps the different libraries in the header and merges the bam files
# The array reads1 have only the fastq from same sample  # TODO: explain better what is 'reads1'
task RunBwaAlignment {

  input {
    Array[String] sampleName
    Array[File] reads
    Array[String] libraries
    ReferenceFasta references
    Int max_cores
    String rm_dupli
  }

  Int disk_size = ceil(size(reads, "GiB") * 2 + size(references.ref_fasta, "GiB") + 20)
  Int memory_size = 4000 * max_cores

  command <<<
    mkdir tmp

    reads_list=( ~{sep=" " reads} )
    lib_list=( ~{sep=" " libraries} )
    sampleName_list=( ~{sep=" " sampleName})
    BAMS=()
    for index in ${!reads_list[*]}; do
      echo "${reads_list[$index]} is in ${lib_list[$index]}"
      bwa_header="@RG\tID:${sampleName_list[$index]}.${lib_list[$index]}\tLB:lib-${lib_list[$index]}\tPL:illumina\tSM:${sampleName_list[$index]}\tPU:FLOWCELL1.LANE1.${lib_list[$index]}"
      /usr/gitc/./bwa mem -t ~{max_cores} -R "${bwa_header}" ~{references.ref_fasta} "${reads_list[$index]}" | \
          java -jar /usr/gitc/picard.jar SortSam \
            I=/dev/stdin \
            O="${sampleName_list[$index]}.${lib_list[$index]}.sorted.bam" \
            TMP_DIR=./tmp \
            SORT_ORDER=coordinate \
            CREATE_INDEX=true;
      mv "${sampleName_list[$index]}.${lib_list[$index]}.sorted.bai" "${sampleName_list[$index]}.${lib_list[$index]}.sorted.bam.bai";
      BAMS+=("I=${sampleName_list[$index]}.${lib_list[$index]}.sorted.bam")
    done

    sampleName_unique=($(echo "${sampleName_list[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

    # Check if there are replicated samples
    for index in ${!sampleName_unique[*]}; do
      NFILES=($(echo ${sampleName_unique[$index]}.*.bam))
      echo ${NFILES[*]}
      echo ${#NFILES[@]}
      REP=()
      if [ "${#NFILES[@]}" -gt 1 ]; then
        for file in ${!NFILES[*]}; do
          REP+=("I=${NFILES[$file]}")
        done
        echo ${REP[*]}

        java -jar /usr/gitc/picard.jar MergeSamFiles ${REP[*]} \
          O=${sampleName_unique[$index]}.sorted_temp.bam \
          CREATE_INDEX=true \
          TMP_DIR=./tmp
      else
        mv ${sampleName_unique[$index]}.*.bam ${sampleName_unique[$index]}.sorted_temp.bam
        mv ${sampleName_unique[$index]}.*.bai ${sampleName_unique[$index]}.sorted_temp.bai
      fi

      if [ "~{rm_dupli}" = "true" ]; then
        java -jar /usr/gitc/picard.jar MarkDuplicates \
            I="${sampleName_unique[$index]}.sorted_temp.bam" \
            O="${sampleName_unique[$index]}.sorted.bam" \
            CLEAR_DT="false" \
            METRICS_FILE= "${sampleName_unique[$index]}_dup_metrics.txt" \
            REMOVE_SEQUENCING_DUPLICATES=true \
            CREATE_INDEX=true
      else
        java -jar /usr/gitc/picard.jar MarkDuplicates \
            I="${sampleName_unique[$index]}.sorted_temp.bam" \
            O="${sampleName_unique[$index]}.sorted_temp2.bam" \
            CLEAR_DT="false" \
            METRICS_FILE= "${sampleName_unique[$index]}_dup_metrics.txt"

        mv "${sampleName_unique[$index]}.sorted_temp.bam" "${sampleName_unique[$index]}.sorted.merged.bam"
        mv "${sampleName_unique[$index]}.sorted_temp.bai" "${sampleName_unique[$index]}.sorted.merged.bai"
      fi

      # Filter by MapQ
      samtools view -bq 10  "${sampleName_unique[$index]}.sorted.merged.bam" >  "${sampleName_unique[$index]}.sorted.merged.filtered.bam"
      samtools index "${sampleName_unique[$index]}.sorted.merged.filtered.bam"
      mv "${sampleName_unique[$index]}.sorted.merged.filtered.bam.bai" "${sampleName_unique[$index]}.sorted.merged.filtered.bai"

    done
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    cpu: max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    preemptible: 3
    # Slurm
    job_name: "RunBwaAlignment"
    mem:"~{memory_size}M"
    time:"05:00:00"
  }

  meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Run [BWA](http://bio-bwa.sourceforge.net/) MEM alignment."
  }

  output {
    Array[File] bam = glob("*.sorted.merged.filtered.bam")
    Array[File] bai = glob("*.sorted.merged.filtered.bai")
    Array[File] dup_metrics = glob("*_dup_metrics.txt")
  }
}

task RunBwaAlignmentSimu {

  input {
    File reads
    Array[File] fastqs
    ReferenceFasta references
    Int max_cores
    String rm_dupli
  }

  Int disk_size = ceil(size(reads, "GiB") + size(fastqs, "GiB") * 2 + size(references.ref_fasta, "GiB"))
  Int memory_size = 14000

  command <<<
    mkdir tmp

    ln -s ~{sep = " " fastqs} .

    for file in $(cat ~{reads}); do

      sample=`basename -s .1.fq $file`
      file_name=`basename $file`

      bwa_header="@RG\tID:${sample}.1\tLB:lib-1\tPL:illumina\tSM:${sample}\tPU:FLOWCELL1.LANE1.1"

      /usr/gitc/./bwa mem -t ~{max_cores} -R "${bwa_header}" ~{references.ref_fasta} "$file_name" | \
          java -jar /usr/gitc/picard.jar SortSam \
            I=/dev/stdin \
            O="${sample}.sorted_temp.bam" \
            TMP_DIR=./tmp \
            SORT_ORDER=coordinate \
            CREATE_INDEX=true

      if [ "~{rm_dupli}" = "true" ]; then
        java -jar /usr/gitc/picard.jar MarkDuplicates \
            I="${sample}.sorted_temp.bam" \
            O="${sample}.sorted.bam" \
            CLEAR_DT="false" \
            METRICS_FILE= "${sample}_dup_metrics.txt" \
            REMOVE_SEQUENCING_DUPLICATES=true \
            CREATE_INDEX=true

      else
        java -jar /usr/gitc/picard.jar MarkDuplicates \
            I="${sample}.sorted_temp.bam" \
            O="${sample}.sorted_temp2.bam" \
            CLEAR_DT="false" \
            METRICS_FILE= "${sample}_dup_metrics.txt"

        mv "${sample}.sorted_temp.bam" "${sample}.sorted.bam"
        mv "${sample}.sorted_temp.bai" "${sample}.sorted.bai"

      fi

      # Filter by MapQ
      # samtools view -bq 10 "${sample}.sorted.bam" > "${sample}.sorted.filtered.bam"
      # samtools index "${sample}.sorted.filtered.bam"
      # mv "${sample}.sorted.filtered.bam.bai" "${sample}.sorted.filtered.bai"

    done

    mkdir dup_metrics
    mv *_dup_metrics.txt dup_metrics
    tar -czvf dup_metrics.tar.gz dup_metrics

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    cpu: max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    preemptible: 3
    # Slurm
    job_name: "RunBwaAlignmentSimu"
    mem:"~{memory_size}M"
    time:"10:00:00"
  }

  meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Run [BWA](http://bio-bwa.sourceforge.net/) MEM alignment in simulated reads."
  }

  output {
    Array[File] bam = glob("*.sorted.bam")
    Array[File] bai = glob("*.sorted.bai")
    File dup_metrics = "dup_metrics.tar.gz"
  }
}

task CreateChunksFastq {
  input {
    Array[String] sampleFile
    Int chunk_size
  }

  Int disk_size = ceil(size(sampleFile, "GiB") * 2)
  Int memory_size = 1000

  command <<<
    set -e
    for i in ~{sep=" " sampleFile}; do echo $i >> lof_sample.txt; done

    split -l ~{chunk_size} lof_sample.txt sample.
  >>>

  runtime {
    docker: "ubuntu:20.04"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "CreateChunksFastq"
    mem:"~{memory_size}M"
    time:"00:05:00"
  }

  meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Split the empirical fastq files into chunks to be aligned in parallel in the next task."
  }

  output {
    Array[File] sample_chunks = glob("sample.*")
  }
}
