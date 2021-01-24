version 1.0

import "snpcalling_empS.wdl"
import "reference_struct.wdl"
import "split_filt_vcf.wdl" as norm_filt
import "CollectAllelicCounts.wdl" as counts


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
    mem:"50GB"
    time:"72:00:00"
    cpu:20
  }

  output {
    File vcf = "freebayes.vcf"
  }
}
