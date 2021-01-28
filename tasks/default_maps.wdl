version 1.0

import "./utilsR.wdl" as utilsR

workflow DefaultMaps {
    input {
     File simu_onemap_obj
     File simu_vcfR
     File vcfR_obj
     File onemap_obj
     File ref_alt_alleles
     File simulated_phases
     String SNPCall_program
     String CountsFrom
     File? multi_obj  
     String seed
     String depth
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
              ref_alt_alleles = ref_alt_alleles,
              simu_onemap_obj = simu_onemap_obj,
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = item.left,
              CountsFrom = CountsFrom,
              simulated_phases = simulated_phases
            }

            call utilsR.ErrorsReport{
              input:
                onemap_obj = select_onemap_obj,
                simu_onemap_obj = simu_onemap_obj,
                SNPCall_program = SNPCall_program,
                GenotypeCall_program = item.left,
                CountsFrom = CountsFrom,
                simu_vcfR = simu_vcfR,
                vcfR_obj = vcfR_obj,
                seed             = seed,
                depth            = depth
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
