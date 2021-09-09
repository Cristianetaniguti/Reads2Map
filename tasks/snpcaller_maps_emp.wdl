version 1.0

import "./utilsR.wdl" as utilsR
import "./utils.wdl" as utils

workflow SNPCallerMaps{
  input {
     File vcf_file
     File? merged_bam
     File? reference
     String cross
     String SNPCall_program
     String GenotypeCall_program
     String CountsFrom
     String parent1
     String parent2
     String chromosome
     String multiallelics
     Int max_cores
    }

  if(multiallelics == 'TRUE') {
      call utils.RunWhatsHap{
          input:
              vcf_file = vcf_file,
              merged_bam = merged_bam,
              reference = reference
      }
  }

  File updated_vcf = select_first([RunWhatsHap.phased_vcf, vcf_file])

  call utilsR.SetProbs{
    input:
      vcf_file = updated_vcf,
      cross = cross,
      parent1 = parent1,
      parent2 = parent2,
      addDefault = TRUE
  }

  Array[String] methods                         = [GenotypeCall_program, GenotypeCall_program + "0.05", GenotypeCall_program + default]
  Array[File] objects                           = [SetProbs.probs_obj_out, SetProbs.globalerror_onemap_obj, SetProbs.default_onemap_obj]
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
        name = "snpcaller_maps"
   }

   output {
      File tar_gz_report = Compress.tar_gz_report
   }

  output{
    File RDatas = MapsReportEmp.maps_RData
    File maps_report = MapsReportEmp.maps_report
    File times = MapsReportEmp.times
    File filters_report = FiltersReportEmp.filters_report
    File errors_report = CheckDepths.errors_report
  }
}

task SNPCallerProbs{
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

      cross <- "~{cross}"

      if(cross == "F1"){
         f1 = NULL
      } else if (cross == "F2"){
         f1 = "F1"
      }

      vcf <- read.vcfR("~{vcf_file}")
      save(vcf, file = "vcfR.RData")

      onemap_obj <- load("~{onemap_obj}")
      onemap_obj <- get(onemap_obj)

      if(any(grepl("freeBayes", vcf@meta))) par <- "GL" else par <- "PL"

      probs <- extract_depth(vcfR.object=vcf,
                               onemap.object=onemap_obj,
                               vcf.par=par,
                               parent1="~{parent1}",
                               parent2="~{parent2}",
                               f1 = f1,
                               recovering=FALSE)

      probs_onemap_obj <- create_probs(onemap.obj = onemap_obj, genotypes_probs=probs)
      save(probs_onemap_obj, file="probs_onemap_obj.RData")

    RSCRIPT

  >>>
  runtime{
    docker:"cristaniguti/reads2map:0.0.1"
    time:"10:00:00"
    mem:"30GB"
    cpu:1
  }

  output{
    File probs_onemap_obj = "probs_onemap_obj.RData"
    File vcfR_obj = "vcfR.RData"
  }
}
