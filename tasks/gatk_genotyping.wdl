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

  call HaplotypeCallerERC {
     input:
       bams = bams,
       bams_index = bais,
       reference_fasta = references.ref_fasta,
       reference_fai = references.ref_fasta_index,
       reference_dict = references.ref_dict
   }

  call norm_filt.SplitFiltVCF {
    input:
      vcf_in=HaplotypeCallerERC.vcf,
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

## Process all samples because it RAD experiments
## usually do not have large ammount of reads (confirm?)
task HaplotypeCallerERC {
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
