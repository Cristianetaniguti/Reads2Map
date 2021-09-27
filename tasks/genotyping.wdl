version 1.0

import "./utilsR.wdl" as utilsR
import "./utils.wdl" as utils

workflow onemapMaps {

  input {
    File vcf_file
    File? merged_bam
    File? reference
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    String cross
    String parent1
    String parent2
    String chromosome
    Int max_cores
    String multiallelics
  }

  call utils.SplitMarkers{
      input:
          vcf_file = vcf_file
  }
  
  call utilsR.ReGenotyping {
      input:
          vcf_file = SplitMarkers.biallelics,
          GenotypeCall_program = GenotypeCall_program,
          cross = cross,
          parent1 = parent1,
          parent2 = parent2,
          max_cores = max_cores
  }

  if (multiallelics == "TRUE") {
    call utils.JointMarkers{
      input:
        biallelic_vcf = ReGenotyping.regeno_vcf,
        multiallelic_vcf = SplitMarkers.multiallelics
    }
  }

  File updated_vcf = select_first([JointMarkers.merged_vcf, ReGenotyping.regeno_vcf])

  call utilsR.SetProbs{
    input:
      vcf_file = updated_vcf,
      cross = cross,
      parent1 = parent1,
      parent2 = parent2,
      multiallelics = multiallelics
  }

  Array[String] methods                         = [GenotypeCall_program, GenotypeCall_program + "0.05"]
  Array[File] objects                           = [SetProbs.probs_onemap_obj, SetProbs.globalerror_onemap_obj]
  Array[Pair[String, File]] methods_and_objects = zip(methods, objects)

  scatter (item in methods_and_objects) {
       call utilsR.CheckDepths {
           input:
              onemap_obj = item.right,
              vcfR_obj = SetProbs.vcfR_obj,
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
        name = "regeno_maps"
   }

   output {
      File tar_gz_report = Compress.tar_gz_report
   }
}
