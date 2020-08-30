version 1.0

import "./utilsR.wdl" as utilsR

workflow DefaultMaps {
    input {
     File simu_onemap_obj
     File onemap_obj
     File tot_mks
     File real_phases
     String SNPCall_program
     String CountsFrom
     String cMbyMb
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
    
         if (defined(multi_obj)) {
            call utilsR.AddMultiallelics{
              input:
                onemap_obj_multi = multi_obj,
                onemap_obj_bi = item.right
           }
         }
        
         File select_onemap_obj = select_first([AddMultiallelics.onemap_obj_both, item.right])        
    
         call utilsR.FiltersReport{
              input:
                onemap_obj = select_onemap_obj,
                SNPCall_program = SNPCall_program,
                GenotypeCall_program = item.left,
                CountsFrom = CountsFrom
          }

          call utilsR.MapsReport{
            input:
              onemap_obj = FiltersReport.onemap_obj_filtered,
              tot_mks = tot_mks,
              simu_onemap_obj = simu_onemap_obj,
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = item.left,
              CountsFrom = CountsFrom,
              cMbyMb = cMbyMb,
              real_phases = real_phases
            }

            call utilsR.ErrorsReport{
              input:
                onemap_obj = select_onemap_obj,
                simu_onemap_obj = simu_onemap_obj,
                SNPCall_program = SNPCall_program,
                GenotypeCall_program = item.left,
                CountsFrom = CountsFrom
            }
     }

     output{
        Array[File] RDatas = MapsReport.maps_RData
        Array[File] maps_report = MapsReport.maps_report
        Array[File] times = MapsReport.times
        Array[File] filters_report = FiltersReport.filters_report
        Array[File] errors_report = ErrorsReport.errors_report

     }

}
