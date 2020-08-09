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
    }

  call GQProbs{
    input:
      vcf_file = vcf_file,
      onemap_obj = onemap_obj,
      cross = cross,
      parent1 = parent1,
      parent2 = parent2
  }

  call utilsR.CheckDepths{
    input:
      onemap_obj = GQProbs.gq_onemap_obj,
      vcfR_obj = GQProbs.vcfR_obj,
      parent1 = parent1,
      parent2 = parent2,
      SNPCall_program = SNPCall_program,
      GenotypeCall_program = GenotypeCall_program,
      CountsFrom = CountsFrom
  }

  call utilsR.FiltersReportEmp{
    input:
      onemap_obj = GQProbs.gq_onemap_obj,
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
      CountsFrom = CountsFrom
  }

  output{
    File RDatas = MapsReportEmp.maps_RData
    File maps_report = MapsReportEmp.maps_report
    File times = MapsReportEmp.times
    File filters_report = FiltersReportEmp.filters_report
    File errors_report = CheckDepths.errors_report
  }
}

task GQProbs{
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

      # MAPS REPORT - GQ
      gq <- extract_depth(vcfR.object=vcf,
                               onemap.object=onemap_obj,
                               vcf.par="GQ",
                               parent1="~{parent1}",
                               parent2="~{parent2}",
                               f1 = f1,
                               recovering=FALSE)

      gq_onemap_obj <- create_probs(onemap.obj = onemap_obj, genotypes_errors=gq)
      save(gq_onemap_obj, file="gq_onemap_obj.RData")

    RSCRIPT

  >>>
  runtime{
    docker:"cristaniguti/onemap_workflows"
    time:"72:00:00"
    # mem:"--nodes=1"
    cpu:1
  }

  output{
    File gq_onemap_obj = "gq_onemap_obj.RData"
    File vcfR_obj = "vcfR.RData"
  }
}
