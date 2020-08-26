version 1.0

import "./utilsR.wdl" as utilsR

workflow DefaultMaps {
    input {
     File onemap_obj
     File vcfR_obj
     String SNPCall_program
     String CountsFrom
     String parent1
     String parent2
     String chromosome
     File? multi_obj
    }

    call utilsR.GlobalError{
      input:
        onemap_obj = onemap_obj
    }

    Array[String] methods                         = ["default", "default0.05"]
    Array[File] objects                           = [onemap_obj, GlobalError.error_onemap_obj]
    Array[Pair[String, File]] methods_and_objects = zip(methods, objects)

    scatter(item in methods_and_objects){

          call utilsR.CheckDepths{
            input:
              onemap_obj = item.right,
              vcfR_obj = vcfR_obj,
              parent1 = parent1,
              parent2 = parent2,
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = item.left,
              CountsFrom = CountsFrom
         }
        
         if (defined(multi_obj)) {
            call utilsR.AddMultiallelics{
              input:
                onemap_obj_multi = multi_obj,
                onemap_obj_bi = item.right
           }
         }
        
         File select_onemap_obj = select_first([AddMultiallelics.onemap_obj_both, item.right])
        
         call utilsR.FiltersReportEmp{
              input:
                onemap_obj = select_onemap_obj,
                SNPCall_program = SNPCall_program,
                GenotypeCall_program = item.left,
                CountsFrom = CountsFrom,
                chromosome = chromosome
          }

          call utilsR.MapsReportEmp{
            input:
              sequence_obj = FiltersReportEmp.onemap_obj_filtered,
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = item.left,
              CountsFrom = CountsFrom
            }
     }

     output{
        Array[File] RDatas = MapsReportEmp.maps_RData
        Array[File] maps_report = MapsReportEmp.maps_report
        Array[File] times = MapsReportEmp.times
        Array[File] filters_report = FiltersReportEmp.filters_report
        Array[File] errors_report = CheckDepths.errors_report
     }
}
