version 1.0

# import "../structs/alignment_struct.wdl"
# import "../structs/reads_simuS.wdl"
import "../structs/snpcalling_empS.wdl"
import "../structs/reference_struct.wdl"
import "./utils.wdl" as utils
import "./utilsR.wdl" as utilsR
import "split_filt_vcf.wdl" as norm_filt


workflow FreebayesGenotyping {
  input {
    Array[File] bam
    Array[File] bai
    Reference references
    String parent1
    String parent2
    String chrom
    String program
    Array[String] sampleNames
    Int max_cores
  }

  call RunFreebayes {
    input:
      reference=references.ref_fasta,
      reference_idx=references.ref_fasta_index,
      bam=bam,
      bai=bai,
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

  call utils.BamCounts {
    input:
      program=program,
      bam=bam,
      ref=references.ref_fasta,
      ref_fai=references.ref_fasta_index,
      ref_dict=references.ref_dict,
      vcf=SplitFiltVCF.vcf_bi,
      tbi=SplitFiltVCF.vcf_bi_tbi
  }


  call utils.BamCounts4Onemap {
    input:
      sampleName=sampleNames,
      counts=BamCounts.counts,
      method = program
  }

  call utilsR.BamDepths2Vcf{
    input:
      vcf_file = SplitFiltVCF.vcf_bi,
      ref_bam = BamCounts4Onemap.ref_bam,
      alt_bam = BamCounts4Onemap.alt_bam,
      example_alleles = BamCounts4Onemap.ref_alt_alleles,
      program = program
  }

  output {
    File vcf_bi = SplitFiltVCF.vcf_bi
    File tbi_bi = SplitFiltVCF.vcf_bi_tbi
    File vcf_multi = SplitFiltVCF.vcf_multi
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
    Int max_cores
  }

  command <<<
   # needed for some singularity versions
   export PATH="/freebayes/vcflib/bin:${PATH}"
   export PATH="/freebayes/scripts:${PATH}"
   export PATH="/freebayes/vcflib/scripts:${PATH}"

   ln -sf ~{sep=" " bam} .
   ln -sf ~{sep=" " bai} .

   freebayes-parallel <(fasta_generate_regions.py ~{reference_idx} 100000) ~{max_cores} \
   --genotype-qualities -f ~{reference}  *.bam > "freebayes.vcf"

  >>>

  runtime {
    docker: "taniguti/freebayes"
    mem:"60GB"
    time:"14:00:00"
    cpu:20
  }

  output {
    File vcf = "freebayes.vcf"
  }
}
