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
    Int chunk_size
  }

  call CreateChunks {
    input:
      bams=bams,
      bams_index=bais,
      chunk_size=chunk_size,
      reference_fasta = references.ref_fasta
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

  Array[String] calling_intervals = read_lines(CreateChunks.interval_list)

  scatter (interval in calling_intervals) {
    call ImportGVCFs {
      input:
        vcfs=flatten(HaplotypeCaller.vcfs),
        vcfs_index=flatten(HaplotypeCaller.vcfs_index),
        reference_fasta=references.ref_fasta,
        reference_fai=references.ref_fasta_index,
        reference_dict=references.ref_dict,
        interval = interval
    }

    call GenotypeGVCFs {
      input:
        workspace_tar = ImportGVCFs.output_workspace,
        interval = interval,
        reference_fasta = references.ref_fasta,
        reference_fai = references.ref_fasta_index,
        reference_dict = references.ref_dict
    }
  }

  call MergeVCFs {
    input:
      input_vcfs = GenotypeGVCFs.vcf,
      input_vcf_indices = GenotypeGVCFs.vcf_tbi,
      ref_fasta = references.ref_fasta,
      ref_fasta_index = references.ref_fasta_index,
      ref_dict = references.ref_dict
  }

  # Simulations
  if(defined(seed)){
    call hard_filt.HardFiltering {
      input:
        references = references,
        vcf_file = MergeVCFs.output_vcf,
        vcf_tbi  = MergeVCFs.output_vcf_index,
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
        vcf_file = MergeVCFs.output_vcf,
        vcf_tbi  = MergeVCFs.output_vcf_index,
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
    File reference_fasta
    Int chunk_size
  }

  command <<<
    set -e
    for i in ~{sep=" " bams}; do echo $i >> lof_bams.txt; done
    for i in ~{sep=" " bams_index}; do echo $i >> lof_bais.txt; done

    split -l ~{chunk_size} lof_bams.txt bams.
    split -l ~{chunk_size} lof_bais.txt bais.

    cat ~{reference_fasta} | grep '>' | tr '\n' ',' | sed '$ s/.$//' | sed 's/,/ \n/g' | sed 's/>//g' > intervals.txt  
  >>>

  runtime {
    docker: "ubuntu:20.04"
    # memory: "2 GB"
    # preemptible: 3
    # cpu: 1
    job_name: "CreateChunks"
    node:"--nodes=1"
    mem:"--mem=1G"
    cpu:"--ntasks=1"
    time:"00:05:00"
  }

  output {
    Array[File] bams_chunks = glob("bams.*")
    Array[File] bais_chunks = glob("bais.*")
    File interval_list = "intervals.txt"
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
      /usr/gitc/gatk4/./gatk --java-options "-Xmx10G -Xms2G" HaplotypeCaller \
        -ERC GVCF \
        -R ~{reference_fasta} \
        -I "$bam" \
        -O "vcfs/${out_name}.g.vcf.gz" \
        --max-reads-per-alignment-start 0 

    done
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    # memory: "4 GB"
    # cpu: 1
    # preemptible: 3
    # disks: "local-disk " + disk_size + " HDD"
    job_name: "HaplotypeCaller"
    node:"--nodes=1"
    mem:"--mem=20G"
    tasks:"--ntasks=1"
    time:"5:00:00"
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

  Int disk_size = ceil(size(vcfs, "GB") + size(reference_fasta, "GB") + 5) * 2

  command <<<
    set -euo pipefail
    grep ">" ~{reference_fasta} | sed 's/^.//' > interval.list
    mkdir gvcfs
    for i in ~{sep=" " vcfs}; do ln -s $i gvcfs/; done
    for i in ~{sep=" " vcfs_index}; do ln -s $i gvcfs/; done

    /usr/gitc/gatk4/./gatk --java-options "-Xmx10G -Xms2G" GenomicsDBImport \
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
    # memory: "4 GB"
    # cpu: 1
    # preemptible: 3
    # disks: "local-disk " + disk_size + " HDD"
    job_name: "ImportGVCFs"
    node:"--nodes=1"
    mem:"--mem=20G"
    tasks:"--ntasks-per-node=5"
    time:"10:00:00"
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

  Int disk_size = ceil(size(reference_fasta, "GB") + 5) * 2

  command <<<
    set -euo pipefail

    tar -xf ~{workspace_tar}

    /usr/gitc/gatk4/./gatk --java-options "-Xmx10G -Xms2G" GenotypeGVCFs \
      -R ~{reference_fasta} \
      -V gendb://cohort_db \
      -L ~{interval} \
      -G StandardAnnotation \
      -O gatk.vcf.gz 

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    # memory: "4 GB"
    # cpu: 1
    # preemptible: 3
    # disks: "local-disk " + disk_size + " HDD"
    job_name: "GenotypeGVCFs"
    node:"--nodes=1"
    mem:"--mem=20G"
    tasks:"--ntasks=1"
    time:"10:00:00"
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
    File ref_fasta
    File ref_fasta_index
    File ref_dict
  }

  command <<<

    /usr/gitc/gatk4/./gatk --java-options "-Xmx10G -Xms2G" \
      MergeVcfs \
      -I ~{sep=' -I' input_vcfs} \
      -O gatk_joint.vcf.gz
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    # memory: "4 GB"
    # cpu: 1
    # preemptible: 3
    # disks: "local-disk " + disk_size + " HDD"
    job_name: "MergeVCFs"
    node:"--nodes=1"
    mem:"--mem=20G"
    tasks:"--ntasks=1"
    time:"10:00:00"
  } 

  output {
    File output_vcf = "gatk_joint.vcf.gz"
    File output_vcf_index = "gatk_joint.vcf.gz.tbi"
  }
}