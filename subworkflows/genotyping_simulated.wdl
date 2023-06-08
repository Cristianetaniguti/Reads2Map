version 1.0

import "../tasks/utilsR.wdl" as utilsR
import "../tasks/utils.wdl" as utils

workflow onemapMaps {

  input {
    File simu_onemap_obj
    File simu_vcfR
    File vcf_file
    File ref_alt_alleles
    File simulated_phases
    String SNPCall_program
    String genotyping_program
    String CountsFrom
    String cross
    Int max_cores
    Int seed
    Int depth
    String multiallelics
    File? multiallelics_file
    Int ploidy
    Float prob_thres
    Array[String] global_errors
    String genoprob_error
    Array[String] genoprob_global_errors
  }


  call utilsR.ReGenotyping {
      input:
          vcf_file = vcf_file,
          GenotypeCall_program = genotyping_program,
          SNPCall_program = SNPCall_program,
          CountsFrom = CountsFrom,
          cross = cross,
          parent1 = "P1",
          parent2 = "P2",
          max_cores = max_cores,
          ploidy = ploidy
  }

  if (multiallelics == "TRUE") {
    call utils.JointMarkers {
      input:
        biallelic_vcf = ReGenotyping.regeno_vcf,
        multiallelic_vcf = multiallelics_file,
        SNPCall_program = SNPCall_program,
        CountsFrom = CountsFrom,
        GenotypeCall_program = genotyping_program
    }
  }

  File updated_vcf = select_first([JointMarkers.merged_vcf, ReGenotyping.regeno_vcf])

  call utilsR.SetProbs {
    input:
      vcf_file = updated_vcf,
      cross = cross,
      parent1 = "P1",
      parent2 = "P2",
      multiallelics = multiallelics,
      SNPCall_program = SNPCall_program,
      global_errors = global_errors,
      genoprob_error = genoprob_error,
      prob_thres = prob_thres,
      genoprob_global_errors = genoprob_global_errors,
      GenotypeCall_program = genotyping_program
  }

  scatter (item in range(length(SetProbs.probs_onemap_obj))) {

       call utilsR.FiltersReport {
            input:
              onemap_obj = SetProbs.probs_onemap_obj[item],
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = SetProbs.probs_onemap_obj_names[item],
              CountsFrom = CountsFrom,
              seed = seed,
              depth = depth
        }

        call utilsR.MapsReport {
          input:
            onemap_obj = SetProbs.probs_onemap_obj[item],
            ref_alt_alleles = ref_alt_alleles,
            simu_onemap_obj = simu_onemap_obj,
            SNPCall_program = SNPCall_program,
            GenotypeCall_program = SetProbs.probs_onemap_obj_names[item],
            CountsFrom = CountsFrom,
            simulated_phases = simulated_phases,
            seed = seed,
            depth = depth,
            max_cores = max_cores
        }

        call utilsR.ErrorsReport {
            input:
              onemap_obj = SetProbs.probs_onemap_obj[item],
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = SetProbs.probs_onemap_obj_names[item],
              CountsFrom = CountsFrom,
              simu_vcfR = simu_vcfR,
              vcfR_obj = SetProbs.vcfR_obj,
              seed = seed,
              depth = depth,
              max_cores = max_cores
        }

   }

   call utils.Compress {
      input:
        RDatas = MapsReport.maps_RData,
        maps_report = MapsReport.maps_report,
        times =  MapsReport.times,
        filters_report = FiltersReport.filters_report,
        errors_report = ErrorsReport.errors_report,
        name = "regeno_maps"
   }

   output {
      File tar_gz_report = Compress.tar_gz_report
   }
}
