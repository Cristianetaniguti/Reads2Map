version 1.0

import "./utilsR.wdl" as utilsR

workflow SupermassaMaps {
  input{
    File simu_onemap_obj
    File onemap_obj
    File vcf_file
    String SNPCall_program
    String CountsFrom
    String cross
    File real_phases
    File tot_mks
    String cMbyMb
  }

  call SupermassaProbs {
    input:
      vcf_file = vcf_file,
      onemap_obj = onemap_obj,
      cross=cross
  }

  call utilsR.GlobalError {
    input:
      onemap_obj = SupermassaProbs.supermassa_onemap_obj
  }

  Array[String] methods                         = ["supermassa", "supermassa0.05"]
  Array[File] objects                           = [SupermassaProbs.supermassa_onemap_obj, GlobalError.error_onemap_obj]
  Array[Pair[String, File]] methods_and_objects = zip(methods, objects)

  scatter(item in methods_and_objects){
       call utilsR.FiltersReport{
            input:
              onemap_obj = item.right,
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
              onemap_obj = item.right,
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = item.left,
              CountsFrom = CountsFrom,
              simu_onemap_obj = simu_onemap_obj
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

task SupermassaProbs{
  input{
    File vcf_file
    File onemap_obj
    String cross
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

       onemap_obj_temp <- load("~{onemap_obj}")
       onemap_obj <- get(onemap_obj_temp)

       supermassa_onemap_obj <- supermassa_genotype(vcfR.object=vcf,
                                    onemap.object=onemap_obj,
                                    vcf.par="AD",
                                    parent1="P1",
                                    parent2="P2",
                                    f1 = f1,
                                    recovering=TRUE,
                                    mean_phred=20,
                                    cores=20,
                                    depths=NULL,
                                    global_error = NULL,
                                    use_genotypes_errors = FALSE,
                                    use_genotypes_probs = TRUE)

       save(supermassa_onemap_obj, file="supermassa_onemap_obj.RData")

     RSCRIPT

  >>>

  runtime{
    docker:"gcr.io/taniguti-backups/onemap:v1"
    time:"72:00:00"
    mem:"--nodes=1"
    cpu:20
  }

  output{
    File supermassa_onemap_obj = "supermassa_onemap_obj.RData"
  }
}
