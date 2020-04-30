version 1.0

import "./utilsR.wdl" as utilsR

workflow DefaultMaps {
    input {
     File onemap_obj
     File vcfR_obj
     String SNPCall_program
     String GenotypeCall_program
     String CountsFrom
     String cross
     String parent1
     String parent2
    }
    
    call utilsR.GlobalError{
      input:
        onemap_obj = onemap_obj,
        SNPCall_program = SNPCall_program,
        GenotypeCall_program = GenotypeCall_program,
        CountsFrom = CountsFrom
    }
    
    Array[String] methods                         = ["default", "default0.05"]
    Array[File] objects                           = [onemap_obj, GlobalError.error_onemap_obj]
    Array[Pair[String, File]] methods_and_objects = zip(methods, objects)
    
    scatter(objects in methods_and_objects){
    
          call utilsR.CheckDepths{
            input:
              onemap_obj = objects.right, 
              vcfR_obj = vcfR_obj, 
              parent1 = parent1,
              parent2 = parent2,
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = objects.left,
              CountsFrom = CountsFrom
        }
         call utilsR.FiltersReport{
              input:
                onemap_obj = objects.right,
                SNPCall_program = SNPCall_program,
                GenotypeCall_program = objects.left,
                CountsFrom = CountsFrom,
                which_workflow = "empirical"
          }
            
          call utilsR.MapsReportEmp{
            input:
              sequence_obj = FiltersReport.onemap_obj_filtered,
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = objects.left,
              CountsFrom = CountsFrom
            }
     }
     
     output{
        Array[File] RDatas = MapsReportEmp.maps_RData
        Array[File] maps_report = MapsReportEmp.maps_report
        Array[File] times = MapsReportEmp.times
        Array[File] filters_report = FiltersReport.filters_report
        Array[File] errors_report = CheckDepths.errors_report
     }
}
