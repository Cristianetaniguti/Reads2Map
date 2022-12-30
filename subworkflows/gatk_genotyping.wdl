version 1.0

import "../structs/dna_seq_structs.wdl"
import "../tasks/chunk_lists.wdl"
import "../tasks/gatk.wdl"
import "../tasks/utils.wdl" as utils

import "../subworkflows/norm_filt_vcf.wdl" as norm_filt
import "../subworkflows/hard_filtering_simulated.wdl" as hard_filt
import "../subworkflows/hard_filtering_empirical.wdl" as hard_filt_emp
import "../subworkflows/MCHap.wdl" as MCHapWf


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
    Boolean mchap
    Boolean replaceAD
    Boolean hardfilters
    Int max_cores
    File? merged_bams
    String? P1
    String? P2
  }

  call chunk_lists.CreateChunksBam {
    input:
      bams=bams,
      bams_index=bais,
      chunk_size=chunk_size,
      reference_fasta = references.ref_fasta
  }

  scatter (chunk in zip(CreateChunksBam.bams_chunks, CreateChunksBam.bais_chunks)) {

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

  Array[String] calling_intervals = read_lines(CreateChunksBam.interval_list)

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
      input_vcf_indices = GenotypeGVCFs.vcf_tbi
  }

  if(hardfilters){
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
  
    File sele_vcf = select_first([HardFiltering.filt_vcf, HardFilteringEmp.filt_vcf])
    File QualPlots = select_first([HardFiltering.Plots, HardFilteringEmp.Plots])

  }

  File filt_vcf = select_first([sele_vcf, MergeVCFs.output_vcf])

  call norm_filt.Normalization {
    input:
      vcf_in= filt_vcf,
      vcf_simu = vcf_simu,
      reference = references.ref_fasta,
      reference_idx = references.ref_fasta_index,
      reference_dict = references.ref_dict,
      program = program,
      counts_source = "vcf",
      ploidy = ploidy
  }

  Map[String, Array[File]] map_bams = {"bam": bams, "bai": bais}

  if(replaceAD) {
    call utils.ReplaceAD {
      input:
        ref_fasta = references.ref_fasta,
        ref_index = references.ref_fasta_index,
        bams = map_bams["bam"],
        bais = map_bams["bai"],
        vcf = Normalization.vcf_norm,
        tbi = Normalization.vcf_norm_tbi,
        program = program,
        counts_source = "bam"
    }
  }

 # MCHap: micro-haplotyping
  if(mchap) {
    call MCHapWf.MCHap {
      input:
          reference = references.ref_fasta,
          reference_idx = references.ref_fasta_index,
          vcf_file = Normalization.vcf_norm,
          n_nodes = 10,
          max_cores = max_cores,
          bams = map_bams["bam"],
          bais = map_bams["bai"],
          ploidy = ploidy,
          merged_bams = merged_bams,
          P1 = P1,
          P2 = P2
    }  

   File vcf_norm_mchap = MCHap.haplo_vcf_merged
 }

 Array[File] gatk_vcfs = select_all([Normalization.vcf_norm, ReplaceAD.bam_vcf]) 
 Array[String] gatk_software = select_all([Normalization.software, ReplaceAD.software])
 Array[String] gatk_counts_source = select_all([Normalization.source, ReplaceAD.source])

  output {
    Array[File] vcfs = gatk_vcfs
    Array[String] vcfs_software = gatk_software
    Array[String] vcfs_counts_source = gatk_counts_source
    File? vcf_multi = vcf_norm_mchap
    File vcfEval = Normalization.vcfEval
    File? Plots = QualPlots
  }
}
