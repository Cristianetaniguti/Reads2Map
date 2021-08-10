version 1.0

import "./utilsR.wdl" as utilsR

workflow SNPCallerMaps{
  input {
     File onemap_obj
     File vcf_file
     String cross
     String SNPCall_program
     String GenotypeCall_program
     String CountsFrom
     String parent1
     String parent2
     String chromosome
     File? multi_obj
     String multiallelics
     Int max_cores
    }

  call SNPCallerProbs{
    input:
      vcf_file = vcf_file,
      onemap_obj = onemap_obj,
      cross = cross,
      parent1 = parent1,
      parent2 = parent2
  }

  call utilsR.CheckDepths{
    input:
      onemap_obj = SNPCallerProbs.probs_onemap_obj,
      vcfR_obj = SNPCallerProbs.vcfR_obj,
      parent1 = parent1,
      parent2 = parent2,
      SNPCall_program = SNPCall_program,
      GenotypeCall_program = GenotypeCall_program,
      CountsFrom = CountsFrom,
      max_cores = max_cores
  }

  if (multiallelics == "TRUE") {
     call utilsR.AddMultiallelics{
         input:
           onemap_obj_multi = multi_obj,
           onemap_obj_bi = SNPCallerProbs.probs_onemap_obj
      }
  }
        
  File select_onemap_obj = select_first([AddMultiallelics.onemap_obj_both, SNPCallerProbs.probs_onemap_obj])

  call utilsR.FiltersReportEmp{
    input:
      onemap_obj = select_onemap_obj,
      SNPCall_program = SNPCall_program,
      GenotypeCall_program = "SNPCaller",
      CountsFrom = "vcf",
      chromosome = chromosome
  }

  call utilsR.MapsReportEmp{
    input:
      sequence_obj = FiltersReportEmp.onemap_obj_filtered,
      SNPCall_program = SNPCall_program,
      GenotypeCall_program = "SNPCaller",
      CountsFrom = CountsFrom,
      max_cores = max_cores
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
