version 1.0

import "./utilsR.wdl" as utilsR

workflow SupermassaMaps{
  input{
    File onemap_obj
    File vcf_file
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    String cross
    String parent1
    String parent2
    String chromosome
  }

  call SupermassaProbs{
    input:
      vcf_file = vcf_file,
      onemap_obj = onemap_obj,
      cross=cross,
      parent1 = parent1,
      parent2 = parent2
  }
  
  call utilsR.GlobalError{
    input:
      onemap_obj = SupermassaProbs.supermassa_onemap_obj,
      SNPCall_program = SNPCall_program,
      GenotypeCall_program = GenotypeCall_program,
      CountsFrom = CountsFrom
  }

  Array[String] methods                         = ["supermassa", "supermassa0.05"]
  Array[File] objects                           = [SupermassaProbs.supermassa_onemap_obj, GlobalError.error_onemap_obj]
  Array[Pair[String, File]] methods_and_objects = zip(methods, objects)
    
  scatter(objects in methods_and_objects){
       call utilsR.CheckDepths{
           input:
              onemap_obj = objects.right, 
              vcfR_obj = SupermassaProbs.vcfR_obj, 
              parent1 = parent1,
              parent2 = parent2,
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = objects.left,
              CountsFrom = CountsFrom
       }
       
       call utilsR.FiltersReportEmp{
            input:
              onemap_obj = objects.right,
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = objects.left,
              CountsFrom = CountsFrom,
              chromosome = chromosome
        }
            
        call utilsR.MapsReportEmp{
          input:
            sequence_obj = FiltersReportEmp.onemap_obj_filtered,
            SNPCall_program = SNPCall_program,
            GenotypeCall_program = objects.left,
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

task SupermassaProbs{
  input{
    File vcf_file
    File onemap_obj
    String cross
    String parent1
    String parent2
  }
  
  command <<<
     R --vanilla --no-save <<RSCRIPT
       library(onemap)
       library(vcfR)
       library(genotyping4onemap)
       
       cross <- "~{cross}"
          
       if(cross == "F1"){
          cross <- "outcross"
          f1 = NULL
       } else if (cross == "F2"){
          cross <- "f2 intercross"
          f1 = "F1"
       }
       
       vcf <- read.vcfR("~{vcf_file}")
       save(vcf, file = "vcfR.RData")
       
       onemap_obj_temp <- load("~{onemap_obj}")
       onemap_obj <- get(onemap_obj_temp)
       
       supermassa_onemap_obj <- supermassa_genotype(vcfR.object=vcf,
                                    onemap.object=onemap_obj,
                                    vcf.par="AD",
                                    parent1="~{parent1}",
                                    parent2="~{parent2}",
                                    f1 = f1,
                                    recovering=TRUE,
                                    mean_phred=20,
                                    cores=3,
                                    depths=NULL,
                                    global_error = NULL,
                                    use_genotypes_errors = FALSE,
                                    use_genotypes_probs = TRUE)
       
       save(supermassa_onemap_obj, file="supermassa_onemap_obj.RData")
  
     RSCRIPT
  >>>
  
  runtime{
    docker:"taniguti/onemap"
    time:"120:00:08"
    mem:"--nodes=1"
    cpu:20
  }
  
  output{
    File supermassa_onemap_obj = "supermassa_onemap_obj.RData"
    File vcfR_obj = "vcfR.RData"
  }
}
