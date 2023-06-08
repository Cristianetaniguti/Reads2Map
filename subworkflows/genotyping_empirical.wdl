version 1.0

import "../tasks/utilsR.wdl" as utilsR
import "../tasks/utils.wdl" as utils

workflow onemapMapsEmp {

  input {
    File vcf_file
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    String cross
    String parent1
    String parent2
    String chromosome
    Int max_cores
    Boolean multiallelics
    File? multiallelics_file
    Int ploidy
    Float prob_thres 
    Array[String] global_errors
    Boolean genoprob_error
    Array[String] genoprob_global_errors
  }

  call utilsR.ReGenotyping {
      input:
          vcf_file = vcf_file,
          SNPCall_program = SNPCall_program,
          CountsFrom = CountsFrom,
          GenotypeCall_program = GenotypeCall_program,
          cross = cross,
          parent1 = parent1,
          parent2 = parent2,
          max_cores = max_cores,
          ploidy = ploidy
  }

  call utils.ApplyRandomFiltersArray as FilterBi{
            input:
                vcfs = [ReGenotyping.regeno_vcf],
                vcfs_SNPCall_software = [SNPCall_program],
                vcfs_Counts_source = [CountsFrom],
                vcfs_GenoCall_software = [GenotypeCall_program + "_biallelic"],
                chromosome = chromosome
  }

  if (multiallelics) {

    Array[File] array_vcf = select_all([multiallelics_file])

    call utils.ApplyRandomFiltersArray as FilterMulti{
            input:
                vcfs = array_vcf,
                vcfs_SNPCall_software = [SNPCall_program],
                vcfs_Counts_source = [CountsFrom],
                vcfs_GenoCall_software = [GenotypeCall_program + "_multiallelic"],
                chromosome = chromosome
    }

    call utils.JointMarkers {
      input:
        biallelic_vcf = FilterBi.vcfs_filt[0],
        multiallelic_vcf = FilterMulti.vcfs_filt[0],
        SNPCall_program = SNPCall_program,
        CountsFrom = CountsFrom,
        GenotypeCall_program = GenotypeCall_program
    }
  }

  File updated_vcf = select_first([JointMarkers.merged_vcf, FilterBi.vcfs_filt[0]])

  call utilsR.SetProbs {
    input:
      vcf_file = updated_vcf,
      cross = cross,
      parent1 = parent1,
      parent2 = parent2,
      multiallelics = multiallelics,
      SNPCall_program = SNPCall_program,
      global_errors = global_errors,
      genoprob_error = genoprob_error,
      prob_thres = prob_thres,
      genoprob_global_errors = genoprob_global_errors,
      GenotypeCall_program = GenotypeCall_program
  }

  scatter (item in range(length(SetProbs.probs_onemap_obj))) {
       call utilsR.CheckDepths {
           input:
              onemap_obj = SetProbs.probs_onemap_obj[item],
              vcfR_obj = SetProbs.vcfR_obj,
              parent1 = parent1,
              parent2 = parent2,
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = SetProbs.probs_onemap_obj_names[item],
              CountsFrom = CountsFrom,
              max_cores = max_cores
       }

       call utilsR.FiltersReportEmp {
            input:
              onemap_obj = SetProbs.probs_onemap_obj[item],
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = SetProbs.probs_onemap_obj_names[item],
              CountsFrom = CountsFrom,
              chromosome = chromosome
        }

        call utilsR.MapsReportEmp {
          input:
            sequence_obj = FiltersReportEmp.onemap_obj_filtered,
            SNPCall_program = SNPCall_program,
            GenotypeCall_program = SetProbs.probs_onemap_obj_names[item],
            CountsFrom = CountsFrom,
            max_cores = max_cores
          }
   }

   call utils.Compress {
      input:
        RDatas = MapsReportEmp.maps_RData,
        maps_report = MapsReportEmp.maps_report,
        times = MapsReportEmp.times,
        filters_report = FiltersReportEmp.filters_report,
        errors_report = CheckDepths.errors_report,
        name = "regeno_maps"
   }

   output {
      File tar_gz_report = Compress.tar_gz_report
      Array[File] maps_report = MapsReportEmp.maps_report
      File regeno_vcf = ReGenotyping.regeno_vcf
   }
}
