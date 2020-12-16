version 1.0

import "../structs/alignment_struct.wdl"
import "../structs/reference_struct.wdl"

# This task considers that is it possible to receive more than one fastq file by sample
# It keeps the different libraries in the header and merges the bam files
# The array reads1 have only the fastq from same sample
task RunBwaAlignment {

  input {
    String sampleName
    Array[File] reads1
    Array[String] libraries
    Reference references
  }

  command <<<
    mkdir tmp
    export PATH=$PATH:/bin
    export PATH=$PATH:/picard.jar

    reads_list=( ~{sep=" " reads1} )
    lib_list=( ~{sep=" " libraries} )
    BAMS=()
    for index in ${!reads_list[*]}; do
      echo "${reads_list[$index]} is in ${lib_list[$index]}"
      bwa_header="@RG\tID:~{sampleName}.${lib_list[$index]}\tLB:lib-${lib_list[$index]}\tPL:illumina\tSM:~{sampleName}\tPU:FLOWCELL1.LANE1.${lib_list[$index]}"
      bwa mem -t 20 -R "${bwa_header}" ~{references.ref_fasta} "${reads_list[$index]}" | \
          java -jar /picard.jar SortSam \
            I=/dev/stdin \
            O="~{sampleName}.${lib_list[$index]}.sorted.bam" \
            TMP_DIR=./tmp \
            SORT_ORDER=coordinate \
            CREATE_INDEX=true;
      mv "~{sampleName}.${lib_list[$index]}.sorted.bai" "~{sampleName}.${lib_list[$index]}.sorted.bam.bai";
      BAMS+=("I=~{sampleName}.${lib_list[$index]}.sorted.bam")
    done

    if [ "${#BAMS[@]}" -gt 1 ]; then
      java -jar /picard.jar MergeSamFiles ${BAMS[@]} \
        O=~{sampleName}.sorted.bam \
        CREATE_INDEX=true \
        TMP_DIR=./tmp
      mv ~{sampleName}.sorted.bai ~{sampleName}.sorted.bam.bai
    else
      mv ~{sampleName}*.bam ~{sampleName}.sorted.bam
      mv ~{sampleName}*.bai ~{sampleName}.sorted.bam.bai
    fi

  >>>

  runtime {
    docker: "kfdrc/bwa-picard:latest-dev"
    time:"72:00:00"
    mem:"--nodes=1"
    cpu:20
  }

  output {
    Alignment algn = {"bam": "~{sampleName}.sorted.bam", "bai": "~{sampleName}.sorted.bam.bai", "sample": "~{sampleName}"}
    File bam = "~{sampleName}.sorted.bam"
    File bai = "~{sampleName}.sorted.bam.bai"
  }
}


# reads1 receive only one fastq
task RunBwaAlignmentSimu {

  input {
    File read1
    Reference references
  }

  command <<<
    mkdir tmp
    export PATH=$PATH:/bin
    export PATH=$PATH:/picard.jar

    filename=`basename '~{read1}'`
    echo $filename
    sample=$(echo $filename | awk -F '[.]' '{print $1}')
    echo $sample > sample.txt
    echo $sample

    bwa_header="@RG\tID:$sample.1\tLB:lib-1\tPL:illumina\tSM:$sample\tPU:FLOWCELL1.LANE1.1"
    bwa mem -t 20 -R "${bwa_header}" ~{references.ref_fasta} ~{read1} | \
        java -jar /picard.jar SortSam \
          I=/dev/stdin \
          O="file.sorted.bam" \
          TMP_DIR=./tmp \
          SORT_ORDER=coordinate \
          CREATE_INDEX=true;
    mv "file.sorted.bai" "file.sorted.bam.bai";

    
  >>>

  runtime {
    docker: "kfdrc/bwa-picard:latest-dev"
    time:"72:00:00"
    mem:"--nodes=1"
    cpu:20
  }

  output {
    Alignment algn = {"bam": "file.sorted.bam", "bai": "file.sorted.bam.bai", "sample": read_string("sample.txt")}
    File bam = "file.sorted.bam"
    File bai = "file.sorted.bam.bai"
  }
}


# Add info to alignment header
task AddAlignmentHeader {
  input {
    String sampleName
    File bam_file
    File bam_idx
  }

  command <<<
    mkdir tmp
    java -jar /gatk/picard.jar AddOrReplaceReadGroups \
      I=~{bam_file} \
      O=~{sampleName}_rg.bam \
      RGLB=lib-~{sampleName} \
      RGPL=illumina \
      RGID=FLOWCELL1.LANE1.~{sampleName} \
      RGSM=~{sampleName} \
      RGPU=FLOWCELL1.LANE1.~{sampleName} \
      CREATE_INDEX=true \
      TMP_DIR=tmp

    mv ~{sampleName}_rg.bai ~{sampleName}_rg.bam.bai

  >>>

  runtime {
    docker: "taniguti/gatk-picard"
    time:"24:00:00"
    mem:"--nodes=1"
    cpu:1
  }

  output {
    Alignment algn = {"bam": "${sampleName}_rg.bam", "bai": "${sampleName}_rg.bam.bai", "sample": "${sampleName}"}
    File bam = "~{sampleName}_rg.bam"
    File bai = "~{sampleName}_rg.bam.bai"
  }
}
