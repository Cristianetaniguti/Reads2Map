version 1.0

import "../../structs/dna_seq_structs.wdl"

import "../../subworkflows/create_alignment_from_families_files.wdl" as fam
import "../../subworkflows/gatk_genotyping.wdl" as gatk
import "../../subworkflows/freebayes_genotyping.wdl" as freebayes


workflow SNPCalling {

  input {
    File samples_info
    ReferenceFasta references
    Int max_cores
    Int chunk_size
    Boolean rm_dupli
    String? P1
    String? P2
    Boolean gatk_mchap
    Boolean hardfilters
    Boolean replaceAD
    Boolean run_gatk
    Boolean run_freebayes
    Int ploidy
    Int n_chrom
  }

  call fam.CreateAlignmentFromFamilies {
    input:
      families_info=samples_info,
      references=references,
      max_cores = max_cores,
      rm_dupli = rm_dupli,
      chunk_size = chunk_size,
      gatk_mchap = gatk_mchap
  }

  if(run_gatk){
    call gatk.GatkGenotyping {
      input:
        bams=CreateAlignmentFromFamilies.bam,
        bais=CreateAlignmentFromFamilies.bai,
        references=references,
        chunk_size = chunk_size,
        ploidy = ploidy,
        program="gatk",
        max_cores = max_cores,
        merged_bams = CreateAlignmentFromFamilies.merged_bam,
        P1 = P1,
        P2 = P2,
        mchap = gatk_mchap,
        hardfilters = hardfilters,
        replaceAD = replaceAD
    }
  }

  if(run_freebayes){
    call freebayes.FreebayesGenotyping {
      input:
        merged_bam=CreateAlignmentFromFamilies.merged_bam,
        references=references,
        program="freebayes",
        max_cores = max_cores,
        ploidy = ploidy,
        replaceAD = replaceAD,
        n_chrom = n_chrom
    }
  }

  Array[File] vcfs_sele = flatten(select_all([GatkGenotyping.vcfs, FreebayesGenotyping.vcfs]))
  Array[File] software_sele = flatten(select_all([GatkGenotyping.vcfs_software, FreebayesGenotyping.vcfs_software]))
  Array[File] source_sele = flatten(select_all([GatkGenotyping.vcfs_counts_source, FreebayesGenotyping.vcfs_counts_source]))

  output {
    Array[File] vcfs = vcfs_sele
    Array[File] vcfs_software = software_sele
    Array[File] vcfs_counts_source = source_sele
    File? gatk_multi_vcf = GatkGenotyping.vcf_multi
    File? gatk_vcfEval = GatkGenotyping.vcfEval
    File? Plots = GatkGenotyping.Plots
    File? freebayes_vcfEval = FreebayesGenotyping.vcfEval
    File? merged_bam = CreateAlignmentFromFamilies.merged_bam
  }
}
