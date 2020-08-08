version 1.0

import "../structs/alignment_struct.wdl"
import "../structs/reads_simuS.wdl"
import "../structs/snpcalling_empS.wdl"
import "./utils.wdl" as utils
import "./utilsR.wdl" as utilsR
import "split_filt_vcf.wdl" as norm_filt


workflow FreebayesGenotyping {
  input {
    Array[Alignment] alignments
    Array[File] bam
    Array[File] bai
    ReferenceFasta references
    SplitVCF splitvcf
    String program
    Array[String] sampleNames
  }

  call RunFreebayes {
    input:
      reference=references.ref_fasta,
      reference_idx=references.ref_fasta_index,
      bam=bam,
      bai=bai
  }

  call norm_filt.SplitFiltVCF{
    input:
      vcf_in=RunFreebayes.vcf,
      program=program,
      chromosome = splitvcf.chromosome,
      reference = references.ref_fasta,
      reference_idx = references.ref_fasta_index,
      parent1 = splitvcf.parent1,
      parent2 = splitvcf.parent2
  }

  scatter (alignment in alignments) {
    call utils.BamCounts {
      input:
        sample=alignment.sample,
        program=program,
        bam=alignment.bam,
        bai=alignment.bai,
        ref=references.ref_fasta,
        ref_fai=references.ref_fasta_index,
        ref_dict=references.ref_dict,
        vcf=SplitFiltVCF.vcf_bi_chr_norm,
        tbi=SplitFiltVCF.vcf_bi_chr_norm_tbi
    }
  }
  
  call utils.BamCounts4Onemap {
    input:
      sampleName=sampleNames,
      counts=BamCounts.counts,
      method = program
  }

  call utilsR.BamDepths2Vcf{
    input:
      vcf_file = SplitFiltVCF.vcf_bi_chr_norm,
      ref_bam = BamCounts4Onemap.ref_bam,
      alt_bam = BamCounts4Onemap.alt_bam,
      example_alleles = BamCounts4Onemap.ref_alt_alleles,
      program = program
  }

  output {
    File vcf = SplitFiltVCF.vcf_bi_chr_norm
    File tbi = SplitFiltVCF.vcf_bi_chr_norm_tbi
    File vcf_bi_tot = SplitFiltVCF.vcf_bi_norm
    File vcf_multi_tot = SplitFiltVCF.vcf_multi_norm
    File vcf_bi_bam_counts = BamDepths2Vcf.bam_vcf
    File alt_bam = BamCounts4Onemap.alt_bam
    File ref_bam = BamCounts4Onemap.ref_bam
  }
}

task RunFreebayes {

  input {
    File reference
    File reference_idx
    Array[File] bam
    Array[File] bai
  }

  command <<<
   freebayes-parallel <(fasta_generate_regions.py ~{reference_idx} 100000) 20 \
   --genotype-qualities -f ~{reference}  ~{sep=" " bam} > "freebayes.vcf"
   
  >>>

  runtime {
    docker: "taniguti/freebayes"
    mem:"--nodes=1"
    time:"72:00:00"
    cpu:20
  }

  output {
    File vcf = "freebayes.vcf"
  }
}
