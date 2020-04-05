version 1.0

import "./utilsR.wdl" as utilsR

workflow UpdogMaps{
  input{
    File simu_onemap_obj
    File onemap_obj
    File vcf_file
    File tot_mks
    File real_phases
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    String cMbyMb
    String cross
  }
  
  call UpdogProbs{
    input:
      vcf_file = vcf_file,
      onemap_obj = onemap_obj,
      cross = cross
  }
  
  call utilsR.GlobalError{
    input:
      onemap_obj = UpdogProbs.updog_onemap_obj,
      SNPCall_program = SNPCall_program,
      GenotypeCall_program = GenotypeCall_program,
      CountsFrom = CountsFrom
  }

  Array[String] methods                         = ["updog", "updog0.05"]
  Array[File] objects                           = [UpdogProbs.updog_onemap_obj, GlobalError.error_onemap_obj]
  Array[Pair[String, File]] methods_and_objects = zip(methods, objects)
    
  scatter(objects in methods_and_objects){
       call utilsR.FiltersReport{
            input:
              onemap_obj = objects.right,
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = objects.left,
              CountsFrom = CountsFrom
        }
            
        call utilsR.MapsReport{
          input:
            onemap_obj = FiltersReport.onemap_obj_filtered,
            tot_mks = tot_mks,
            simu_onemap_obj = simu_onemap_obj,
            SNPCall_program = SNPCall_program,
            GenotypeCall_program = objects.left,
            CountsFrom = CountsFrom,
            cMbyMb = cMbyMb,
            real_phases = real_phases
          }
            
          call utilsR.ErrorsReport{
            input:
              onemap_obj = objects.right,
              simu_onemap_obj = simu_onemap_obj,
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = objects.left,
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

task UpdogProbs{
  input{
    File vcf_file
    File onemap_obj
    String cross
  }
  
  command <<<
     R --vanilla --no-save <<RSCRIPT
       library(onemap)
       library(vcfR)
       
       vcf <- read.vcfR("~{vcf_file}")
      
       cross <- "~{cross}"
          
       if(cross == "F1"){
          cross <- "outcross"
          f1 = NULL
       } else if (cross == "F2"){
          cross <- "f2 intercross"
          f1 = "F1"
       }

       
       onemap_obj_temp <- load("~{onemap_obj}")
       onemap_obj <- get(onemap_obj_temp)
       
       updog_onemap_obj <- updog_genotype(vcfR.object=vcf,
                                    onemap.object=onemap_obj,
                                    vcf.par="AD",
                                    parent1="P1",
                                    parent2="P2",
                                    f1 = f1,
                                    recovering=TRUE,
                                    mean_phred=20,
                                    cores=3,
                                    depths=NULL,
                                    global_error = NULL,
                                    use_genotypes_errors = FALSE,
                                    use_genotypes_probs = TRUE)
       
       save(updog_onemap_obj, file="updog_onemap_obj.RData")
  
     RSCRIPT
  >>>
  
  runtime{
    docker:"taniguti/onemap"
  }
  
  output{
    File updog_onemap_obj = "updog_onemap_obj.RData"
  }
}