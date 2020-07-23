version 1.0

import "../structs/alignment_struct.wdl"

task RunBwaAlignment {

  input {
    String sampleName
    File ref
    Array[File] reads1
    Array[String] libraries
    File geno_amb
    File geno_ann
    File geno_bwt
    File geno_pac
    File geno_sa
  }

  command <<<
    echo ~{geno_amb} ~{geno_ann} ~{geno_bwt} ~{geno_pac} ~{geno_sa}
    mkdir tmp
    export PATH=$PATH:/bin
    export PATH=$PATH:/picard.jar

    reads_list=( ~{sep=" " reads1} )
    lib_list=( ~{sep=" " libraries} )
    BAMS=()
    for index in ${!reads_list[*]}; do
      echo "${reads_list[$index]} is in ${lib_list[$index]}"
      bwa_header="@RG\tID:~{sampleName}.${lib_list[$index]}\tLB:lib-${lib_list[$index]}\tPL:illumina\tSM:~{sampleName}\tPU:FLOWCELL1.LANE1.${lib_list[$index]}"
      bwa mem -t 20 -R "${bwa_header}" ~{ref} "${reads_list[$index]}" | \
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
