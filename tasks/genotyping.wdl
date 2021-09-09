version 1.0

import "./utilsR.wdl" as utilsR
import "./utils.wdl" as utils

workflow OneMapMaps {

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

  call utilsR.ReGenotyping {
      input:
          vcf_file = vcf_file,
          GenotypeCall_program = GenotypeCall_program,
          cross = dataset.cross,
          parent1 = dataset.parent1,
          parent2 = dataset.parent2,
          max_cores = max_cores
  
  if(multiallelics == 'TRUE') {
      call utils.RunWhatsHap{
          input:
              vcf_file = ReGenotyping.regeno_vcf,
              merged_bam = merged_bam,
              reference = reference
      }
  }

  File updated_vcf = select_first([RunWhatsHap.phased_vcf, ReGenotyping.regeno_vcf])

  call utilsR.SetProbs{
    input:
      vcf_file = updated_vcf,
      cross = cross,
      parent1 = parent1,
      parent2 = parent2,
      addDefault = FALSE
  }

  Array[String] methods                         = [GenotypeCall_program, GenotypeCall_program + "0.05"]
  Array[File] objects                           = [SetProbs.probs_obj_out, SetProbs.globalerror_onemap_obj]
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

   call Compress {
      input:
        RDatas = MapsReportEmp.maps_RData
        maps_report = MapsReportEmp.maps_report
        times = MapsReportEmp.times
        filters_report = FiltersReportEmp.filters_report
        errors_report = CheckDepths.errors_report,
        name = "regeno_maps"
   }

   output {
      File tar_gz_report = Compress.tar_gz_report
   }
}
