version 1.0

import "../structs/snpcalling_empS.wdl"
import "../structs/reference_struct.wdl"
import "split_filt_vcf.wdl" as norm_filt
import "utils.wdl" as utils
import "hard_filtering.wdl" as hard_filt
import "hard_filtering_emp.wdl" as hard_filt_emp


workflow GatkGenotyping {
  input {
    Array[File] bams
    Array[File] bais
    Reference references
    String program
    File? vcf_simu
    Int? seed
    Int? depth
  }

  call CreateChunks {
    input:
      bams=bams,
      bams_index=bais,
      chunk_size=30
  }

  scatter (chunk in zip(CreateChunks.bams_chunks, CreateChunks.bais_chunks)) {

    call HaplotypeCaller {
      input:
        bams = read_lines(chunk.left),
        bams_index = read_lines(chunk.right),
        reference_fasta = references.ref_fasta,
        reference_fai = references.ref_fasta_index,
        reference_dict = references.ref_dict
    }
  }

  call GATKJointCall {
    input:
      vcfs=flatten(HaplotypeCaller.vcfs),
      vcfs_index=flatten(HaplotypeCaller.vcfs_index),
      reference_fasta=references.ref_fasta,
      reference_fai=references.ref_fasta_index,
      reference_dict=references.ref_dict
  }

  # Simulations
  if(defined(seed)){
    call hard_filt.HardFiltering {
      input:
        references = references,
        vcf_file = GATKJointCall.vcf,
        vcf_tbi  = GATKJointCall.vcf_tbi,
        simu_vcf = vcf_simu,
        seed = seed,
        depth = depth
    }
  }

  # Empirical
  if(!defined(seed)){
    call hard_filt_emp.HardFilteringEmp {
      input:
        references = references,
        vcf_file = GATKJointCall.vcf,
        vcf_tbi  = GATKJointCall.vcf_tbi
    }
  }

  File filt_vcf = select_first([HardFiltering.filt_vcf, HardFilteringEmp.filt_vcf])
  File QualPlots = select_first([HardFiltering.Plots, HardFilteringEmp.Plots])

  call norm_filt.Normalization {
    input:
      vcf_in= filt_vcf,
      vcf_simu = vcf_simu,
      program=program,
      reference = references.ref_fasta,
      reference_idx = references.ref_fasta_index,
      reference_dict = references.ref_dict
  }

  Map[String, Array[File]] map_bams = {"bam": bams, "bai": bais}

  call utils.ReplaceAD {
    input:
      ref_fasta = references.ref_fasta,
      ref_index = references.ref_fasta_index,
      bams = map_bams["bam"],
      bais = map_bams["bai"],
      vcf = Normalization.vcf_norm,
      tbi = Normalization.vcf_norm_tbi,
      program = program
  }

  output {
    File vcf_norm = Normalization.vcf_norm
    File vcf_norm_bamcounts = ReplaceAD.bam_vcf
    File vcfEval = Normalization.vcfEval
    File Plots = QualPlots
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
    # memory: "2 GB"
    # preemptible: 3
    # cpu: 1
    job_name: "CreateChunks"
    node:"--nodes=1"
    mem:"--mem=1GB"
    cpu:"--ntasks=1"
    time:"00:05:00"
  }

  output {
    Array[File] bams_chunks = glob("bams.*")
    Array[File] bais_chunks = glob("bais.*")
  }
}

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
  }

  Int disk_size = ceil(size(reference_fasta, "GB") + size(bams, "GB") * 2)

  command <<<
    set -euo pipefail

    for bam in ~{sep=" " bams}; do ln -s $bam .; done
    for bai in ~{sep=" " bams_index}; do ln -s $bai .; done

    mkdir vcfs
    ## gvcf for each sample
    for bam in *.bam; do
      out_name=$(basename -s ".bam" "$bam")
      /gatk/gatk HaplotypeCaller \
        -ERC GVCF \
        -R ~{reference_fasta} \
        -I "$bam" \
        -O "vcfs/${out_name}.g.vcf.gz" \
        --max-reads-per-alignment-start 0
    done
  >>>

  runtime {
    docker: "taniguti/gatk-picard"
    # memory: "4 GB"
    # cpu: 1
    # preemptible: 3
    # disks: "local-disk " + disk_size + " HDD"
    job_name: "HaplotypeCaller"
    node:"--nodes=1"
    mem:"--mem=32GB"
    cpu:"--ntasks=1"
    time:"05:00:00"
  }

  output {
    Array[File] vcfs = glob("vcfs/*.vcf.gz")
    Array[File] vcfs_index = glob("vcfs/*.vcf.gz.tbi")
  }
}

task GATKJointCall {
  input {
    Array[File] vcfs
    Array[File] vcfs_index
    File reference_fasta
    File reference_fai
    File reference_dict
  }

  Int disk_size = ceil(size(vcfs, "GB") + size(reference_fasta, "GB") + 5) * 2

  command <<<
    set -euo pipefail
    grep ">" ~{reference_fasta} | sed 's/^.//' > interval.list
    mkdir gvcfs
    for i in ~{sep=" " vcfs}; do ln -s $i gvcfs/; done
    for i in ~{sep=" " vcfs_index}; do ln -s $i gvcfs/; done

    gatk --java-options "-Xmx3700m -Xms2g" GenomicsDBImport \
      --batch-size 50 \
      --genomicsdb-workspace-path cohort_db \
      -L interval.list \
      -V $(find gvcfs/*.g.vcf.gz -type l | paste -d',' -s | sed 's/,/ -V /g')

    gatk --java-options "-Xmx3700m -Xms2g" GenotypeGVCFs \
      -R ~{reference_fasta} \
      -V gendb://cohort_db \
      -G StandardAnnotation \
      -O gatk.vcf.gz

  >>>

  runtime {
    docker: "taniguti/gatk-picard"
    # memory: "4 GB"
    # cpu: 1
    # preemptible: 3
    # disks: "local-disk " + disk_size + " HDD"
    job_name: "GATKJointCall"
    node:"--nodes=1"
    mem:"--mem=64GB"
    cpu:"--ntasks=1"
    time:"05:00:00"
  }

  output {
    File vcf = "gatk.vcf.gz"
    File vcf_tbi = "gatk.vcf.gz.tbi"
  }
}
