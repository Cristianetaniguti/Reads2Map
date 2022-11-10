version 1.0

import "../structs/struct_reference.wdl"
import "./norm_filt_vcf.wdl" as norm_filt
import "./utils.wdl" as utils


workflow FreebayesGenotyping {
  input {
    Array[File] bams
    Array[File] bais
    Reference references
    String program
    Int max_cores
    File? vcf_simu
  }

  call RunFreebayes {
    input:
      reference=references.ref_fasta,
      reference_idx=references.ref_fasta_index,
      bam=bams,
      bai=bais,
      max_cores = max_cores
  }

  call norm_filt.Normalization {
    input:
      vcf_in = RunFreebayes.vcf,
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
  }
}

task RunFreebayes {

  input {
    File reference
    File reference_idx
    Array[File] bam
    Array[File] bai
    Int max_cores
  }

  Int disk_size = ceil(size(reference, "GiB") + size(bam, "GiB") +  50)
  Int memory_size = ceil(size(bam, "MiB") * 1.25)

  command <<<
   # needed for some singularity versions
   ln -sf ~{sep=" " bam} .
   ln -sf ~{sep=" " bai} .

   freebayes-parallel <(fasta_generate_regions.py ~{reference_idx} 100000) ~{max_cores} \
   --genotype-qualities -f ~{reference} *.bam > "freebayes.vcf"

  >>>

  runtime {
    docker: "cristaniguti/freebayes:0.0.1"
    cpu: max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "RunFreebayes"
    mem:"~{memory_size}M"
    time:"48:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Split genomic regions and runs [freebayes](https://github.com/freebayes/freebayes) parallelized."
  }

  output {
    File vcf = "freebayes.vcf"
  }
}
