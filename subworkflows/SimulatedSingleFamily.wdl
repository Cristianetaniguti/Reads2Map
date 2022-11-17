version 1.0

import "../structs/population_structs.wdl"

import "../tasks/utils.wdl" as utils
import "../tasks/utilsR.wdl" as utilsR
import "../tasks/JointReports.wdl" as reports

import "../subworkflows/create_alignment_from_read_simulations.wdl" as simulation
import "../subworkflows/genotyping_simulated.wdl" as genotyping
import "../subworkflows/gusmap_maps_simulated.wdl" as gusmap
import "../subworkflows/snpcaller_maps_simulated.wdl" as snpcaller
import "../subworkflows/freebayes_genotyping.wdl" as freebayes
import "../subworkflows/gatk_genotyping.wdl" as gatk


workflow SimulatedSingleFamily {

  input {
    ReferenceFasta references
    Family family
    Sequencing sequencing
    String? filters
    Int max_cores
    Int chunk_size
    Int ploidy
    Boolean gatk_mchap
    Boolean hardfilters
    Boolean replaceAD
    Int n_chrom
  }

  call simulation.CreateAlignmentFromSimulation {
    input:
      sequencing = sequencing,
      references=references,
      family=family,
      max_cores = max_cores,
      chunk_size = chunk_size
  }

  call gatk.GatkGenotyping {
    input:
      bams       = CreateAlignmentFromSimulation.bam,
      bais       = CreateAlignmentFromSimulation.bai,
      references = references,
      program    = "gatk",
      vcf_simu   = CreateAlignmentFromSimulation.true_vcf,
      seed       = family.seed,
      depth      = sequencing.depth,
      chunk_size = chunk_size,
      ploidy     = ploidy,
      mchap      = gatk_mchap,
      max_cores  = max_cores,
      merged_bams = CreateAlignmentFromSimulation.merged_bam,
      P1 = "P1",
      P2 = "P2",
      hardfilters = hardfilters,
      replaceAD = replaceAD
  }

  call freebayes.FreebayesGenotyping {
    input:
      merged_bam   = CreateAlignmentFromSimulation.merged_bam,
      references = references,
      program    = "freebayes",
      max_cores  = max_cores,
      vcf_simu   = CreateAlignmentFromSimulation.true_vcf,
      ploidy     = family.ploidy,
      replaceAD  = replaceAD,
      n_chrom    = n_chrom
  }

  call utilsR.vcf2onemap as truth_vcf {
    input:
      vcf_file = CreateAlignmentFromSimulation.true_vcf,
      cross = family.cross,
      parent1 = "P1",
      parent2 = "P2"
  }

  if (defined(filters)) {
      call utils.ApplyRandomFilters {
          input:
              gatk_vcf = GatkGenotyping.vcf_norm,
              freebayes_vcf = FreebayesGenotyping.vcf_norm,
              gatk_vcf_bam_counts = GatkGenotyping.vcf_norm_bamcounts,
              freebayes_vcf_bam_counts = FreebayesGenotyping.vcf_norm_bamcounts,
              filters = filters,
              chromosome = sequencing.chromosome
      }
  }

  File filtered_gatk_vcf = select_first([ApplyRandomFilters.gatk_vcf_filt,  GatkGenotyping.vcf_norm])
  File filtered_gatk_vcf_bamcounts = select_first([ApplyRandomFilters.gatk_vcf_bam_counts_filt, GatkGenotyping.vcf_norm_bamcounts])
  File filtered_freebayes_vcf = select_first([ApplyRandomFilters.freebayes_vcf_filt, FreebayesGenotyping.vcf_norm])
  File filtered_freebayes_vcf_bamcounts = select_first([ApplyRandomFilters.freebayes_vcf_bam_counts_filt, FreebayesGenotyping.vcf_norm_bamcounts])

  call utils.GetMarkersPos {
    input:
      true_vcf = CreateAlignmentFromSimulation.true_vcf,
      filtered_gatk_vcf = filtered_gatk_vcf,
      filtered_gatk_vcf_bamcounts = filtered_gatk_vcf_bamcounts,
      filtered_freebayes_vcf = filtered_freebayes_vcf,
      filtered_freebayes_vcf_bamcounts = filtered_freebayes_vcf_bamcounts,
      depth = sequencing.depth,
      seed = family.seed
  }

  PopulationAnalysis gatk_processing = {"method": "gatk", "vcf": filtered_gatk_vcf, "bam": filtered_gatk_vcf_bamcounts}
  PopulationAnalysis freebayes_processing = {"method": "freebayes", "vcf": filtered_freebayes_vcf, "bam": filtered_freebayes_vcf_bamcounts}

  scatter (analysis in [gatk_processing, freebayes_processing]){

    Map[String, File] vcfs = {"vcf": analysis.vcf, "bam": analysis.bam}

    scatter (origin in ["vcf", "bam"]){

        call utils.SplitMarkers as splitgeno {
             input:
                vcf_file = vcfs[origin]
        }

        # Suggestion for better SuperMASSA, updog and polyRAD performance
        call utilsR.FilterSegregation {
             input:
                vcf_file = splitgeno.biallelics,
                parent1 = "P1",
                parent2 = "P2"
        }

        call genotyping.onemapMaps as updogMaps {
          input:
            simu_onemap_obj = truth_vcf.onemap_obj,
            vcf_file = FilterSegregation.vcf_filtered,
            #vcf_file = splitgeno.biallelics,
            genotyping_program = "updog",
            ref_alt_alleles = CreateAlignmentFromSimulation.ref_alt_alleles,
            simulated_phases = CreateAlignmentFromSimulation.simulated_phases,
            SNPCall_program = analysis.method,
            CountsFrom = origin,
            cross = family.cross,
            max_cores = max_cores,
            simu_vcfR = truth_vcf.vcfR_obj,
            seed = family.seed,
            depth = sequencing.depth,
            multiallelics = sequencing.multiallelics,
            multiallelics_file = splitgeno.multiallelics
        }

        call genotyping.onemapMaps as supermassaMaps {
          input:
            simu_onemap_obj = truth_vcf.onemap_obj,
            vcf_file = FilterSegregation.vcf_filtered,
            #vcf_file = splitgeno.biallelics,
            genotyping_program = "supermassa",
            ref_alt_alleles = CreateAlignmentFromSimulation.ref_alt_alleles,
            simulated_phases = CreateAlignmentFromSimulation.simulated_phases,
            SNPCall_program = analysis.method,
            CountsFrom = origin,
            cross = family.cross,
            max_cores = max_cores,
            simu_vcfR = truth_vcf.vcfR_obj,
            seed = family.seed,
            depth = sequencing.depth,
            multiallelics = sequencing.multiallelics,
            multiallelics_file = splitgeno.multiallelics
        }

        call genotyping.onemapMaps as polyradMaps {
          input:
            simu_onemap_obj = truth_vcf.onemap_obj,
            vcf_file = FilterSegregation.vcf_filtered,
            #vcf_file = splitgeno.biallelics,
            genotyping_program = "polyrad",
            ref_alt_alleles = CreateAlignmentFromSimulation.ref_alt_alleles,
            simulated_phases = CreateAlignmentFromSimulation.simulated_phases,
            SNPCall_program = analysis.method,
            CountsFrom = origin,
            cross = family.cross,
            max_cores = max_cores,
            simu_vcfR = truth_vcf.vcfR_obj,
            seed = family.seed,
            depth = sequencing.depth,
            multiallelics = sequencing.multiallelics,
            multiallelics_file = splitgeno.multiallelics
        }
    }

    call utils.SplitMarkers as splitvcf {
         input:
           vcf_file = analysis.vcf
    }

    call utils.SplitMarkers as splitbam {
         input:
           vcf_file = analysis.bam
    }

    call gusmap.gusmapMaps {
      input:
        simu_onemap_obj = truth_vcf.onemap_obj,
        vcf_file = splitvcf.biallelics,
        new_vcf_file = splitbam.biallelics,
        SNPCall_program = analysis.method,
        GenotypeCall_program = "gusmap",
        ref_alt_alleles = CreateAlignmentFromSimulation.ref_alt_alleles,
        simulated_phases = CreateAlignmentFromSimulation.simulated_phases,
        seed = family.seed,
        depth = sequencing.depth,
        max_cores = max_cores
    }

    call snpcaller.SNPCallerMaps {
      input:
        simu_onemap_obj = truth_vcf.onemap_obj,
        vcf_file = splitvcf.biallelics,
        ref_alt_alleles = CreateAlignmentFromSimulation.ref_alt_alleles,
        simulated_phases = CreateAlignmentFromSimulation.simulated_phases,
        cross = family.cross,
        SNPCall_program = analysis.method,
        GenotypeCall_program = "SNPCaller",
        CountsFrom = "vcf",
        simu_vcfR = truth_vcf.vcfR_obj,
        seed = family.seed,
        depth = sequencing.depth,
        max_cores = max_cores,
        multiallelics = sequencing.multiallelics,
        multiallelics_file = splitvcf.multiallelics,
        multiallelics_mchap = GatkGenotyping.vcf_multi,
        mchap = gatk_mchap
    }
  }

  # Compress files
  call reports.JointReportsSimu {
    input:
      SNPCaller                 = SNPCallerMaps.tar_gz_report,
      updog                     = flatten(updogMaps.tar_gz_report),
      polyrad                   = flatten(polyradMaps.tar_gz_report),
      supermassa                = flatten(supermassaMaps.tar_gz_report),
      gusmap_files              = gusmapMaps.tar_gz_report,
      GATK_eval                 = GatkGenotyping.vcfEval,
      Freebayes_eval            = FreebayesGenotyping.vcfEval,
      multiallelics_file        = splitvcf.multiallelics,
      max_cores                 = max_cores,
      seed                      = family.seed,
      depth                     = sequencing.depth
  }

  output {
    File data1_depths_geno_prob   = JointReportsSimu.data1_depths_geno_prob
    File data2_maps               = JointReportsSimu.data2_maps
    File data3_filters            = JointReportsSimu.data3_filters
    File data4_times              = JointReportsSimu.data4_times
    File data5_SNPCall_efficiency = JointReportsSimu.data5_SNPCall_efficiency
    File data6_RDatas             = JointReportsSimu.data6_RDatas
    File data7_gusmap             = JointReportsSimu.data7_gusmap
    File data8_names              = JointReportsSimu.data8_names
    File data10_counts            = JointReportsSimu.data10_counts
    File simu_haplo               = CreateAlignmentFromSimulation.simu_haplo
    File? Plots                    = GatkGenotyping.Plots
    File positions                = GetMarkersPos.positions
  }
}
