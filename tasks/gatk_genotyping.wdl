version 1.0

import "../structs/alignment_struct.wdl"
import "../structs/reads_simuS.wdl"
import "./utils.wdl" as utils

workflow GatkGenotyping {
  input {
    Array[Alignment] alignments
    ReferenceFasta references
    String program

  }

  scatter (alignment in alignments) {
    call HaplotypeCallerERC {
      input:
        ref        = references.ref_fasta,
        geno_fai   = references.ref_fasta_index,
        sampleName = alignment.sample,
        bam_rg     = alignment.bam,
        bam_rg_idx = alignment.bai,
        geno_dict  = references.ref_dict
    }

    call GenotypeGVCFs {
      input:
      sample=alignment.sample,
      vcf_in=HaplotypeCallerERC.GVCF,
      fasta=references.ref_fasta,
      fasta_fai=references.ref_fasta_index,
      fasta_dict=references.ref_dict
    }
  }

  call utils.BcftoolsMerge {
    input:
      prefix=program,
      vcfs=GenotypeGVCFs.vcf,
      tbis=GenotypeGVCFs.tbi
  }

  call utils.VcftoolsApplyFilters {
    input:
      vcf_in=BcftoolsMerge.vcf,
      max_missing=0.75,
      min_alleles=2,
      max_alleles=2,
      maf=0.05,
      program=program
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
        vcf=VcftoolsApplyFilters.vcf,
        tbi=VcftoolsApplyFilters.tbi
    }
  }

  output {
    File vcf = VcftoolsApplyFilters.vcf
    File tbi = VcftoolsApplyFilters.tbi
    Array[File] counts = BamCounts.counts
  }
}


task HaplotypeCallerERC {
  input {
    File ref
    File geno_fai
    String sampleName
    File bam_rg
    File bam_rg_idx
    File geno_dict
  }

  command <<<
    /gatk/gatk HaplotypeCaller \
      -ERC GVCF \
      -R ~{ref} \
      -I ~{bam_rg} \
      -O ~{sampleName}_rawLikelihoods.g.vcf \
      --max-reads-per-alignment-start 0
  >>>

  runtime {
    docker: "taniguti/gatk-picard"
  }

  output {
    File GVCF = "${sampleName}_rawLikelihoods.g.vcf"
    File GVCF_idx = "${sampleName}_rawLikelihoods.g.vcf.idx"
  }
}


# Variant calling on gVCF
task GenotypeGVCFs {

  input {
    String sample
    File vcf_in
    File fasta
    File fasta_fai
    File fasta_dict
  }

  command <<<
    /gatk/gatk GenotypeGVCFs \
        -R ~{fasta} \
        -O ~{sample}.gatk.vcf.gz \
        --variant ~{vcf_in}
  >>>

  runtime {
    docker: "taniguti/gatk-picard"
  }

  output {
    File vcf = "~{sample}.gatk.vcf.gz"
    File tbi = "~{sample}.gatk.vcf.gz.tbi"
  }
}
