version 1.0

import "../structs/struct_reference.wdl"
import "norm_filt_vcf.wdl" as norm_filt
import "utils.wdl" as utils
import "MCHap.wdl" as MCHapWf
import "hard_filtering-simulated.wdl" as hard_filt
import "hard_filtering-empirical.wdl" as hard_filt_emp


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
    Int ploidy
    String mchap
    Int max_cores
    File? merged_bams
    String? P1
    String? P2
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
        reference_dict = references.ref_dict,
        ploidy = ploidy,
        chunk_size = chunk_size
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

 # MCHap: micro-haplotyping
 if(mchap == "TRUE") {
   
   Array[File] counts_source = [Normalization.vcf_norm, ReplaceAD.bam_vcf]

   scatter (one_vcf in counts_source){
      call MCHapWf.MCHap{
        input:
          reference = references.ref_fasta,
          reference_idx = references.ref_fasta_index,
          vcf_file = one_vcf, 
          n_nodes = 10,
          max_cores = max_cores,
          bams = map_bams["bam"],
          bais = map_bams["bai"],
          ploidy = ploidy, 
          merged_bams = merged_bams,
          P1 = P1,
          P2 = P2
      }
   }

   File vcf_norm_mchap = MCHap.haplo_vcf_merged[0]
   File vcf_bam_mchap = MCHap.haplo_vcf_merged[1]
 }

  output {
    File? vcf_multi = vcf_norm_mchap
    File? vcf_multi_bamcounts = vcf_bam_mchap
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

  Int disk_size = ceil(size(reference_fasta, "GiB") + 2)

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
    cpu: 1
    # Cloud
    memory:"1000 MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "CreateChunks"
    mem:"1G"
    time:"00:05:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Split the the samples BAM alignment files in chunks."
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
    Int ploidy
    Int chunk_size
  }

  Int disk_size = ceil((size(bams, "GiB") + 30) + size(reference_fasta, "GiB")) + 20
  Int memory_size = ceil(10000 * chunk_size)
  Int max_cores = ceil(chunk_size * 4 + 2)

  command <<<
    set -euo pipefail

    for bam in ~{sep=" " bams}; do ln -s $bam .; done
    for bai in ~{sep=" " bams_index}; do ln -s $bai .; done

    mkdir vcfs
    ## gvcf for each sample
    for bam in *.bam; do
      out_name=$(basename -s ".bam" "$bam")
      /usr/gitc/gatk4/./gatk --java-options "-Xms8000m -Xmx9000m" HaplotypeCaller \
        -ERC GVCF \
        -R ~{reference_fasta} \
        -ploidy ~{ploidy} \
        -I "$bam" \
        -O "vcfs/${out_name}.g.vcf.gz" \
        --max-reads-per-alignment-start 0 &
    done

    wait
    
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

  command <<<
    set -euo pipefail
    grep ">" ~{reference_fasta} | sed 's/^.//' > interval.list
    mkdir gvcfs
    for i in ~{sep=" " vcfs}; do ln -s $i gvcfs/; done
    for i in ~{sep=" " vcfs_index}; do ln -s $i gvcfs/; done

    /usr/gitc/gatk4/./gatk --java-options "-Xms8000m -Xmx25000m" GenomicsDBImport \
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
    memory:"26000 MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "ImportGVCFs"
    mem:"26000M"
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

  command <<<
    set -euo pipefail

    tar -xf ~{workspace_tar}

    /usr/gitc/gatk4/./gatk --java-options "-Xms8000m -Xmx25000m" GenotypeGVCFs \
      -R ~{reference_fasta} \
      -V gendb://cohort_db \
      -L ~{interval} \
      -G StandardAnnotation \
      -O gatk.vcf.gz 

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    cpu: 2
    # Cloud
    memory:"26000 MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "GenotypeGVCFs"
    mem:"26000M"
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
    File ref_fasta
    File ref_fasta_index
    File ref_dict
  }

  Int disk_size = ceil(size(input_vcfs, "GiB") * 2.5) + 10

  command <<<

    /usr/gitc/gatk4/./gatk --java-options "-Xms2000m -Xmx2500m" \
      MergeVcfs \
        -I ~{sep=' -I' input_vcfs} \
        -O gatk_joint.vcf.gz
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    cpu: 1
    # Cloud
    memory:"3000 MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "MergeVCFs"
    mem:"3000M"
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