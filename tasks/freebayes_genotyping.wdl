version 1.0

import "../structs/snpcalling_empS.wdl"
import "../structs/reference_struct.wdl"
import "./split_filt_vcf.wdl" as norm_filt
import "./utils.wdl" as utils


workflow FreebayesGenotyping {
  input {
    Array[File] bams
    Array[File] bais
    Reference references
    String parent1
    String parent2
    String program
    Array[String] sample_names
    Int max_cores
  }

  call RunFreebayes {
    input:
      reference=references.ref_fasta,
      reference_idx=references.ref_fasta_index,
      bam=bams,
      bai=bais,
      max_cores = max_cores
  }

  call norm_filt.SplitFiltVCF {
    input:
      vcf_in=RunFreebayes.vcf,
      program=program,
      reference = references.ref_fasta,
      reference_idx = references.ref_fasta_index,
      parent1 = parent1,
      parent2 = parent2
  }

  Map[String, Array[File]] map_bams = {"bam": bams, "bai": bais}

  call utils.ReplaceAD {
    input:
      ref_fasta = references.ref_fasta,
      ref_index = references.ref_fasta_index,
      bams = map_bams["bam"],
      bais = map_bams["bai"],
      vcf = SplitFiltVCF.vcf_biallelics,
      tbi = SplitFiltVCF.vcf_biallelics_tbi,
      program = program
  }

  output {
    File vcf_biallelics = SplitFiltVCF.vcf_biallelics
    File tbvcf_biallelics_tbii_bi = SplitFiltVCF.vcf_biallelics_tbi
    File vcf_multiallelics = SplitFiltVCF.vcf_multiallelics
    File vcf_biallelics_bamcounts = ReplaceAD.bam_vcf
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

  Int disk_size = ceil(size(reference, "GB") + size(bam, "GB") * 2)

  command <<<
   # needed for some singularity versions
   export PATH="/freebayes/vcflib/bin:${PATH}"
   export PATH="/freebayes/scripts:${PATH}"
   export PATH="/freebayes/vcflib/scripts:${PATH}"

   ln -sf ~{sep=" " bam} .
   ln -sf ~{sep=" " bai} .

   freebayes-parallel <(fasta_generate_regions.py ~{reference_idx} 100000) ~{max_cores} \
   --genotype-qualities -f ~{reference} *.bam > "freebayes.vcf"

  >>>

  runtime {
    docker: "taniguti/freebayes"
    memory: "4 GB"
    preemptible: 3
    cpu: 4
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File vcf = "freebayes.vcf"
  }
}
