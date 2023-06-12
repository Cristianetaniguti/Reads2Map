version 1.0

import "../tasks/utilsR.wdl" as utilsR
import "../tasks/utils.wdl" as utils

workflow SNPCallerMapsEmp {
  input {
     File vcf_file
     String cross
     String SNPCall_program
     String GenotypeCall_program
     String CountsFrom
     String parent1
     String parent2
     String chromosome
     String multiallelics
     File? multiallelics_file
     Int max_cores
     Float prob_thres 
     Array[String] global_errors
     Boolean genoprob_error
     Array[String] genoprob_global_errors

    }

  if (multiallelics == "TRUE") {
    call utils.JointMarkers {
      input:
        biallelic_vcf = vcf_file,
        multiallelic_vcf = multiallelics_file,
        SNPCall_program = SNPCall_program,
        CountsFrom = CountsFrom,
        GenotypeCall_program = GenotypeCall_program
    }
   }

  File updated_vcf = select_first([JointMarkers.merged_vcf, vcf_file])

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
        name = "snpcaller_maps"
   }

   output {
      File tar_gz_report = Compress.tar_gz_report
   }
}
