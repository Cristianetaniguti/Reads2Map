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
     File? multiallelics_mchap
     String mchap
     Int max_cores

    }

  if (multiallelics == "TRUE") {
    call utils.JointMarkers {
      input:
        biallelic_vcf = vcf_file,
        multiallelic_vcf = multiallelics_file
    }
   }

  File updated_vcf = select_first([JointMarkers.merged_vcf, vcf_file])

  call utilsR.SetProbsDefault {
    input:
      vcf_file = updated_vcf,
      cross = cross,
      parent1 = parent1,
      parent2 = parent2,
      SNPCall_program = SNPCall_program,
      multiallelics = multiallelics,
      multiallelics_mchap = multiallelics_mchap,
      mchap = mchap
  }

  Array[String] methods                         = [GenotypeCall_program, GenotypeCall_program + "0.05", GenotypeCall_program + "default"]
  Array[File] objects                           = [SetProbsDefault.probs_onemap_obj, SetProbsDefault.globalerror_onemap_obj, SetProbsDefault.default_onemap_obj]
  Array[Pair[String, File]] methods_and_objects = zip(methods, objects)

  scatter (item in methods_and_objects) {
       call utilsR.CheckDepths {
           input:
              onemap_obj = item.right,
              vcfR_obj = SetProbsDefault.vcfR_obj,
              parent1 = parent1,
              parent2 = parent2,
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = item.left,
              CountsFrom = CountsFrom,
              max_cores = max_cores
       }

       call utilsR.FiltersReportEmp {
            input:
              onemap_obj = item.right,
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = item.left,
              CountsFrom = CountsFrom,
              chromosome = chromosome
        }

        call utilsR.MapsReportEmp {
          input:
            sequence_obj = FiltersReportEmp.onemap_obj_filtered,
            SNPCall_program = SNPCall_program,
            GenotypeCall_program = item.left,
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
