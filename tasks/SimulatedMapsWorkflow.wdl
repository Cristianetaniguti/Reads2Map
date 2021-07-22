version 1.0


import "./create_alignment_from_read_simulations.wdl" as simulation
import "./gatk_genotyping.wdl" as gatk
import "./freebayes_genotyping.wdl" as freebayes
import "./utils.wdl" as utils
import "./utilsR.wdl" as utilsR
import "./default_maps.wdl" as default
import "./snpcaller_maps.wdl" as snpcaller
import "./gusmap_maps.wdl" as gusmap
import "./genotyping-simulated.wdl" as genotyping


struct PopulationAnalysis {
    String method
    File vcf
    File bam
    File multi
}

workflow SimulatedMapsWorkflow {

  input {
    Reference references
    Family family
    Sequencing sequencing
    String? filters
    Int max_cores
  }

  call simulation.CreateAlignmentFromSimulation {
    input:
      sequencing = sequencing,
      references=references,
      family=family,
      max_cores = max_cores
  }

  call gatk.GatkGenotyping {
    input:
      bams=CreateAlignmentFromSimulation.bam,
      bais=CreateAlignmentFromSimulation.bai,
      references=references,
      program="gatk",
      parent1 = "P1",
      parent2 = "P2",
      vcf_simu = CreateAlignmentFromSimulation.true_vcf,
      seed    = family.seed,
      depth   = sequencing.depth
  }

  call freebayes.FreebayesGenotyping {
    input:
      bams=CreateAlignmentFromSimulation.bam,
      bais=CreateAlignmentFromSimulation.bai,
      references=references,
      program="freebayes",
      parent1 = "P1",
      parent2 = "P2",
      max_cores = max_cores,
      vcf_simu = CreateAlignmentFromSimulation.true_vcf
  }

  call utilsR.vcf2onemap as truth_vcf {
    input:
      vcf_file = CreateAlignmentFromSimulation.true_vcf,
      cross = family.cross,
      SNPCall_program = "simu",
      parent1 = "P1",
      parent2 = "P2"
  }

  if (defined(filters)) {
      call utils.ApplyRandomFilters {
          input:
              gatk_vcf = GatkGenotyping.vcf_biallelics,
              freebayes_vcf = FreebayesGenotyping.vcf_biallelics,
              gatk_vcf_bam_counts = GatkGenotyping.vcf_biallelics_bamcounts,
              freebayes_vcf_bam_counts = FreebayesGenotyping.vcf_biallelics_bamcounts,
              filters = filters,
              chromosome = sequencing.chromosome
      }
  }

    File filtered_gatk_vcf = select_first([ApplyRandomFilters.gatk_vcf_filt,  GatkGenotyping.vcf_biallelics])
    File filtered_gatk_vcf_bamcounts = select_first([ApplyRandomFilters.gatk_vcf_bam_counts_filt, GatkGenotyping.vcf_biallelics_bamcounts])
    File filtered_freebayes_vcf = select_first([ApplyRandomFilters.freebayes_vcf_filt, FreebayesGenotyping.vcf_biallelics])
    File filtered_freebayes_vcf_bamcounts = select_first([ApplyRandomFilters.freebayes_vcf_bam_counts_filt, FreebayesGenotyping.vcf_biallelics_bamcounts])

    PopulationAnalysis gatk_processing = {"multi": GatkGenotyping.vcf_multiallelics, "method": "gatk", "vcf": filtered_gatk_vcf, "bam": filtered_gatk_vcf_bamcounts}
    PopulationAnalysis freebayes_processing = {"multi": FreebayesGenotyping.vcf_multiallelics, "method": "freebayes", "vcf": filtered_freebayes_vcf, "bam": filtered_freebayes_vcf_bamcounts}

  scatter (analysis in [gatk_processing, freebayes_processing]){

    call utilsR.vcf2onemap {
      input:
        vcf_file = analysis.vcf,
        cross = family.cross,
        SNPCall_program = analysis.method,
        parent1 = "P1",
        parent2 = "P2"
    }

    call utilsR.MultiVcf2onemap{
      input:
          multi = analysis.multi,
          cross = family.cross,
          SNPCall_program = analysis.method,
          parent1 = "P1",
          parent2 = "P2",
          seed    = family.seed,
          depth   = sequencing.depth,
          multiallelics = sequencing.multiallelics
    }

    call default.DefaultMaps {
      input:
        onemap_obj = vcf2onemap.onemap_obj,
        simu_onemap_obj = truth_vcf.onemap_obj,
        ref_alt_alleles = CreateAlignmentFromSimulation.ref_alt_alleles,
        simulated_phases = CreateAlignmentFromSimulation.simulated_phases,
        SNPCall_program = analysis.method,
        CountsFrom = "vcf",
        multi_obj = MultiVcf2onemap.onemap_obj,
        simu_vcfR = truth_vcf.vcfR_obj,
        vcfR_obj = vcf2onemap.vcfR_obj,
        seed = family.seed,
        depth = sequencing.depth,
        max_cores = max_cores,
        multiallelics = sequencing.multiallelics
    }

    call snpcaller.SNPCallerMaps{
      input:
        simu_onemap_obj = truth_vcf.onemap_obj,
        onemap_obj = vcf2onemap.onemap_obj,
        vcf_file = analysis.vcf,
        ref_alt_alleles = CreateAlignmentFromSimulation.ref_alt_alleles,
        simulated_phases = CreateAlignmentFromSimulation.simulated_phases,
        cross = family.cross,
        SNPCall_program = analysis.method,
        GenotypeCall_program = "SNPCaller",
        CountsFrom = "vcf",
        multi_obj = MultiVcf2onemap.onemap_obj,
        simu_vcfR = truth_vcf.vcfR_obj,
        seed = family.seed,
        depth = sequencing.depth,
        max_cores = max_cores,
        multiallelics = sequencing.multiallelics
    }

    Map[String, File] vcfs = {"vcf": analysis.vcf, "bam": analysis.bam}

    scatter (origin in ["vcf", "bam"]){
        call genotyping.SnpBasedGenotypingSimulatedMaps as UpdogMaps {
          input:
            simu_onemap_obj = truth_vcf.onemap_obj,
            onemap_obj = vcf2onemap.onemap_obj,
            vcf_file = vcfs[origin],
            genotyping_program = "updog",
            ref_alt_alleles = CreateAlignmentFromSimulation.ref_alt_alleles,
            simulated_phases = CreateAlignmentFromSimulation.simulated_phases,
            SNPCall_program = analysis.method,
            CountsFrom = origin,
            cross = family.cross,
            multi_obj = MultiVcf2onemap.onemap_obj,
            max_cores = max_cores,
            simu_vcfR = truth_vcf.vcfR_obj,
            seed = family.seed,
            depth = sequencing.depth,
            multiallelics = sequencing.multiallelics
        }

        call genotyping.SnpBasedGenotypingSimulatedMaps as SupermassaMaps {
          input:
            simu_onemap_obj = truth_vcf.onemap_obj,
            onemap_obj = vcf2onemap.onemap_obj,
            vcf_file = vcfs[origin],
            genotyping_program = "supermassa",
            ref_alt_alleles = CreateAlignmentFromSimulation.ref_alt_alleles,
            simulated_phases = CreateAlignmentFromSimulation.simulated_phases,
            SNPCall_program = analysis.method,
            CountsFrom = origin,
            cross = family.cross,
            multi_obj = MultiVcf2onemap.onemap_obj,
            max_cores = max_cores,
            simu_vcfR = truth_vcf.vcfR_obj,
            seed = family.seed,
            depth = sequencing.depth,
            multiallelics = sequencing.multiallelics
        }

        call genotyping.SnpBasedGenotypingSimulatedMaps as PolyradMaps {
          input:
            simu_onemap_obj = truth_vcf.onemap_obj,
            onemap_obj = vcf2onemap.onemap_obj,
            vcf_file = vcfs[origin],
            genotyping_program = "polyrad",
            ref_alt_alleles = CreateAlignmentFromSimulation.ref_alt_alleles,
            simulated_phases = CreateAlignmentFromSimulation.simulated_phases,
            SNPCall_program = analysis.method,
            CountsFrom = origin,
            cross = family.cross,
            multi_obj = MultiVcf2onemap.onemap_obj,
            max_cores = max_cores,
            simu_vcfR = truth_vcf.vcfR_obj,
            seed = family.seed,
            depth = sequencing.depth,
            multiallelics = sequencing.multiallelics
        }
      }

      call gusmap.GusmapMaps {
        input:
          simu_onemap_obj = truth_vcf.onemap_obj,
          vcf_file = analysis.vcf,
          new_vcf_file = analysis.bam,
          SNPCall_program = analysis.method,
          GenotypeCall_program = "gusmap",
          ref_alt_alleles = CreateAlignmentFromSimulation.ref_alt_alleles,
          simulated_phases = CreateAlignmentFromSimulation.simulated_phases,
          seed = family.seed,
          depth = sequencing.depth,
          max_cores = max_cores
      }
  }

  call utilsR.JointReports {
    input:
    default_RDatas            = flatten(DefaultMaps.RDatas),
    default_maps_report       = flatten(DefaultMaps.maps_report),
    default_filters_report    = flatten(DefaultMaps.filters_report),
    default_errors_report     = flatten(DefaultMaps.errors_report),
    default_times_report      = flatten(DefaultMaps.times),
    SNPCaller_RDatas          = SNPCallerMaps.RDatas,
    multi_names               = MultiVcf2onemap.multi_names,
    SNPCaller_maps_report     = SNPCallerMaps.maps_report,
    SNPCaller_filters_report  = SNPCallerMaps.filters_report,
    SNPCaller_errors_report   = SNPCallerMaps.errors_report,
    SNPCaller_times_report    = SNPCallerMaps.times,
    Updog_RDatas              = flatten(flatten(UpdogMaps.RDatas)),
    Updog_maps_report         = flatten(flatten(UpdogMaps.maps_report)),
    Updog_filters_report      = flatten(flatten(UpdogMaps.filters_report)),
    Updog_errors_report       = flatten(flatten(UpdogMaps.errors_report)),
    Updog_times_report        = flatten(flatten(UpdogMaps.times)),
    Polyrad_RDatas            = flatten(flatten(PolyradMaps.RDatas)),
    Polyrad_maps_report       = flatten(flatten(PolyradMaps.maps_report)),
    Polyrad_filters_report    = flatten(flatten(PolyradMaps.filters_report)),
    Polyrad_errors_report     = flatten(flatten(PolyradMaps.errors_report)),
    Polyrad_times_report      = flatten(flatten(PolyradMaps.times)),
    Supermassa_RDatas         = flatten(flatten(SupermassaMaps.RDatas)),
    Supermassa_maps_report    = flatten(flatten(SupermassaMaps.maps_report)),
    Supermassa_filters_report = flatten(flatten(SupermassaMaps.filters_report)),
    Supermassa_errors_report  = flatten(flatten(SupermassaMaps.errors_report)),
    Supermassa_times_report   = flatten(flatten(SupermassaMaps.times)),
    Gusmap_RDatas             = flatten(GusmapMaps.RDatas),
    Gusmap_maps_report        = flatten(GusmapMaps.maps_report),
    Gusmap_times_report       = flatten(GusmapMaps.times),
    GATK_eval                 = GatkGenotyping.vcfEval,
    Freebayes_eval            = FreebayesGenotyping.vcfEval,
    max_cores                 = max_cores,
    seed                      = family.seed,
    depth                     = sequencing.depth,
  }

  output {
    File data1_depths_geno_prob   = JointReports.data1_depths_geno_prob
    File data2_maps               = JointReports.data2_maps
    File data3_filters            = JointReports.data3_filters
    File data4_times              = JointReports.data4_times
    File data5_SNPCall_efficiency = JointReports.data5_SNPCall_efficiency
    File data6_RDatas             = JointReports.data6_RDatas
    File data7_gusmap             = JointReports.data7_gusmap
    File data8_names              = JointReports.data8_names
    File data10_counts            = JointReports.data10_counts
    File simu_haplo               = CreateAlignmentFromSimulation.simu_haplo
    File Plots                    = GatkGenotyping.Plots
  }
}
