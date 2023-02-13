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
    Boolean rm_dupli = false
    String? P1
    String? P2
    Boolean gatk_mchap = false
    Boolean hardfilters = true
    Boolean replaceAD = false
    Boolean run_gatk = true
    Boolean run_freebayes = true
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

  Array[Array[File]] vcfs_sele = select_all([GatkGenotyping.vcfs, FreebayesGenotyping.vcfs])
  Array[Array[String]] software_sele = select_all([GatkGenotyping.vcfs_software, FreebayesGenotyping.vcfs_software])
  Array[Array[String]] source_sele = select_all([GatkGenotyping.vcfs_counts_source, FreebayesGenotyping.vcfs_counts_source])

  output {
    Array[File] vcfs = flatten(vcfs_sele)
    Array[String] vcfs_software = flatten(software_sele)
    Array[String] vcfs_counts_source = flatten(source_sele)
    File? gatk_multi_vcf = GatkGenotyping.vcf_multi
    File? gatk_vcfEval = GatkGenotyping.vcfEval
    File? Plots = GatkGenotyping.Plots
    File? freebayes_vcfEval = FreebayesGenotyping.vcfEval
    File? merged_bam = CreateAlignmentFromFamilies.merged_bam
  }
}
