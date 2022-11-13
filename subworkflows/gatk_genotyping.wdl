version 1.0

import "../structs/struct_reference.wdl"
import "../tasks/custom/chunk_lists.wdl"
import "../tasks/gatk.wdl"
import "../tasks/utils.wdl" as utils

import "norm_filt_vcf.wdl" as norm_filt
import "hard_filtering_simulated.wdl" as hard_filt
import "hard_filtering_empirical.wdl" as hard_filt_emp
import "MCHap.wdl" as MCHapWf


workflow GatkGenotyping {
  input {
    Array[File] bams
    Array[File] bais
    ReferenceFasta references
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

  call chunk_lists.CreateChunks {
    input:
      bams=bams,
      bams_index=bais,
      chunk_size=chunk_size,
      reference_fasta = references.ref_fasta
  }

  scatter (chunk in zip(CreateChunks.bams_chunks, CreateChunks.bais_chunks)) {

    call gatk.HaplotypeCaller {
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
    call gatk.ImportGVCFs {
      input:
        vcfs=flatten(HaplotypeCaller.vcfs),
        vcfs_index=flatten(HaplotypeCaller.vcfs_index),
        reference_fasta=references.ref_fasta,
        reference_fai=references.ref_fasta_index,
        reference_dict=references.ref_dict,
        interval = interval
    }

    call gatk.GenotypeGVCFs {
      input:
        workspace_tar = ImportGVCFs.output_workspace,
        interval = interval,
        reference_fasta = references.ref_fasta,
        reference_fai = references.ref_fasta_index,
        reference_dict = references.ref_dict
    }
  }

  call gatk.MergeVCFs {
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
