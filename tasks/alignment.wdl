version 1.0

import "../structs/struct_alignment.wdl"
import "../structs/struct_reference.wdl"

# This task considers that is it possible to receive more than one fastq file by sample
# It keeps the different libraries in the header and merges the bam files
# The array reads1 have only the fastq from same sample
task RunBwaAlignment {

  input {
    Array[String] sampleName
    Array[File] reads
    Array[String] libraries
    Reference references
    Int max_cores
    String rm_dupli
  }

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

      if [ "~{rm_dupli}" = "TRUE" ]; then
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
    done
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    # memory: "1 GB"
    # cpu:4
    # preemptible: 3
    # disks: "local-disk " + 10 + " HDD"
    job_name: "RunBwaAlignment"
    node:"--nodes=1"
    mem:"--mem=32GB"
    tasks:"--ntasks-per-node=11"
    time:"00:40:00"
  }

  output {
    Array[File] bam = glob("*.sorted.merged.bam")
    Array[File] bai = glob("*.sorted.merged.bai")
    Array[File] dup_metrics = glob("*_dup_metrics.txt")
  }
}


task RunBwaAlignmentSimu {

  input {
    Array[File] reads
    Reference references
    Int max_cores
    String rm_dupli
  }

  command <<<
    mkdir tmp
    for file in ~{sep= " " reads}; do

      sample=`basename -s .1.fq $file`

      bwa_header="@RG\tID:${sample}.1\tLB:lib-1\tPL:illumina\tSM:${sample}\tPU:FLOWCELL1.LANE1.1"

      /usr/gitc/./bwa mem -t ~{max_cores} -R "${bwa_header}" ~{references.ref_fasta} $file | \
          java -jar /usr/gitc/picard.jar SortSam \
            I=/dev/stdin \
            O="${sample}.sorted_temp.bam" \
            TMP_DIR=./tmp \
            SORT_ORDER=coordinate \
            CREATE_INDEX=true

      if [ "~{rm_dupli}" = "TRUE" ]; then
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
                  
    done

    mkdir dup_metrics
    mv *_dup_metrics.txt dup_metrics
    tar -czvf dup_metrics.tar.gz dup_metrics

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    # memory: "1 GB"
    # cpu:4
    # preemptible: 3
    # disks: "local-disk " + 10 + " HDD"
    job_name: "RunBwaAlignmentSimu"
    node:"--nodes=1"
    mem:"--mem=1GB"
    tasks:"--ntasks-per-node=11"
    time:"00:20:00"
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

  command <<<
    set -e
    for i in ~{sep=" " sampleFile}; do echo $i >> lof_sample.txt; done

    split -l ~{chunk_size} lof_sample.txt sample.
  >>>

  runtime {
    docker: "ubuntu:20.04"
    # memory: "2 GB"
    # preemptible: 3
    # cpu: 1
    job_name: "CreateChunksFastq"
    node:"--nodes=1"
    mem:"--mem=1G"
    cpu:"--ntasks=1"
    time:"00:05:00"
  }

  output {
    Array[File] sample_chunks = glob("sample.*")
  }
}