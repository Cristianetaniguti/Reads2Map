version 1.0

import "snpcalling_empS.wdl"
import "reference_struct.wdl"
import "split_filt_vcf.wdl" as norm_filt
import "CollectAllelicCounts.wdl" as counts


workflow GatkGenotyping {
  input {
    Array[File] bams
    Array[File] bais
    Reference references
    String program
    String parent1
    String parent2
    Array[String] sample_names
  }

  call CreateChunks {
    input:
      bams=bams,
      bams_index=bais,
      chunk_size=50
  }

  scatter (chunk in zip(CreateChunks.bams_chunks, CreateChunks.bais_chunks)) {

    call HaplotypeCallerJointCall {
      input:
        bams = read_lines(chunk.left),
        bams_index = read_lines(chunk.right),
        reference_fasta = references.ref_fasta,
        reference_fai = references.ref_fasta_index,
        reference_dict = references.ref_dict
    }
  }

  call GatherVCFs {
    input:
      input_vcfs=HaplotypeCallerJointCall.vcf,
      input_vcfs_indexes=HaplotypeCallerJointCall.vcf_index,
      output_vcf_name="gatk.vcf.gz"
  }

  call norm_filt.SplitFiltVCF {
    input:
      vcf_in=GatherVCFs.output_vcf,
      program=program,
      reference = references.ref_fasta,
      reference_idx = references.ref_fasta_index,
      parent1 = parent1,
      parent2 = parent2
  }

  call counts.CollectAllelicCountsToVcf {
    input:
      program=program,
      sample_names=sample_names,
      bams=bams,
      bams_index=bais,
      references=references,
      vcf_biallelics_splitted=SplitFiltVCF.vcf_biallelics,
      vcf_biallelics_tbi_splitted=SplitFiltVCF.vcf_biallelics_tbi
  }

  output {
    File vcf_biallelics = SplitFiltVCF.vcf_biallelics
    File vcf_biallelics_tbi = SplitFiltVCF.vcf_biallelics_tbi
    File vcf_multiallelics = SplitFiltVCF.vcf_multiallelics
    File vcf_biallelics_bamcounts = CollectAllelicCountsToVcf.vcf_biallelics_bamcounts
    File alt_bam = CollectAllelicCountsToVcf.alt_bam
    File ref_bam = CollectAllelicCountsToVcf.ref_bam
  }
}

task CreateChunks {
  input {
    Array[String] bams
    Array[String] bams_index
    Int chunk_size
  }

  command <<<
    set -e
    for i in ~{sep=" " bams}; do echo $i >> lof_bams.txt; done
    for i in ~{sep=" " bams_index}; do echo $i >> lof_bais.txt; done

    split -l ~{chunk_size} lof_bams.txt bams.
    split -l ~{chunk_size} lof_bais.txt bais.
  >>>

  runtime {
    docker: "ubuntu:20.04"
  }

  output {
    Array[File] bams_chunks = glob("bams.*")
    Array[File] bais_chunks = glob("bais.*")
  }
}

## Process all samples because it RAD experiments
## usually do not have large ammount of reads (confirm?)
task HaplotypeCallerJointCall {
  input {
    File reference_fasta
    File reference_dict
    File reference_fai
    Array[File] bams
    Array[File] bams_index
  }

  command <<<
    set -euo pipefail

    mkdir vcfs
    ## gvcf for each sample
    for bam in ~{sep= " " bams}; do
      out_name=$(basename $bam)
      /gatk/gatk HaplotypeCaller \
        -ERC GVCF \
        -R ~{reference_fasta} \
        -I "$bam" \
        -O "vcfs/${out_name}.rawLikelihoods.g.vcf.gz" \
        --max-reads-per-alignment-start 0
    done

    grep ">" ~{reference_fasta} | sed 's/^.//' > interval.list
    # Create string as " sample.vcf.gz -V sample2.vcf.gz -V ..."
    vcfs=$(find vcfs/ -maxdepth 1 -name '*.vcf.gz'| awk '{$1=$1}1' OFS=' -V ' RS='')
    ## Combine into genomic database
    /gatk/gatk GenomicsDBImport \
        --genomicsdb-workspace-path gatk_database \
        -L interval.list \
        -V $vcfs

    /gatk/gatk GenotypeGVCFs \
        -R ~{reference_fasta} \
        -O gatk.vcf.gz \
        -G StandardAnnotation \
        -V gendb://gatk_database
  >>>

  runtime {
    docker: "taniguti/gatk-picard"
    mem:"10GB"
    cpu:1
    time:"120:00:00"
  }

  output {
    File vcf = "gatk.vcf.gz"
    File vcf_index = "gatk.vcf.gz.tbi"
  }
}


task GatherVCFs {

  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_vcf_name
  }
  # using MergeVcfs instead of GatherVcfs so we can create indices
  # WARNING	2015-10-28 15:01:48	GatherVcfs	Index creation not currently supported when gathering block compressed VCFs.
  command {
    java -Xmx2g -jar /usr/picard/picard.jar \
      MergeVcfs \
      INPUT=~{sep=' INPUT=' input_vcfs} \
      OUTPUT=~{output_vcf_name}
  }
  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.24.1"
    memory: "3 GB"
    disks: "local-disk " + 10 + " HDD"
    preemptible: 3
  }
}
