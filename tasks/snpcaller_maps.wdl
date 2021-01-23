version 1.0

import "./utilsR.wdl" as utilsR

workflow SNPCallerMaps{
  input {
     File simu_onemap_obj
     File onemap_obj
     File vcf_file
     File ref_alt_alleles
     File simulated_phases
     String cross
     String SNPCall_program
     String GenotypeCall_program
     String CountsFrom
     File? multi_obj
    }


  call GQProbs{
    input:
      vcf_file = vcf_file,
      onemap_obj = onemap_obj,
      cross = cross
  }
  
  if (defined(multi_obj)) {
      call utilsR.AddMultiallelics{
          input:
            onemap_obj_multi = multi_obj,
            onemap_obj_bi = GQProbs.gq_onemap_obj
      }
  }
        
  File select_onemap_obj = select_first([AddMultiallelics.onemap_obj_both, GQProbs.gq_onemap_obj])

  call utilsR.FiltersReport{
    input:
      onemap_obj = select_onemap_obj,
      SNPCall_program = SNPCall_program,
      GenotypeCall_program = GenotypeCall_program,
      CountsFrom = "vcf"
  }

  call utilsR.MapsReport{
    input:
      onemap_obj = FiltersReport.onemap_obj_filtered,
      ref_alt_alleles = ref_alt_alleles,
      simu_onemap_obj = simu_onemap_obj,
      SNPCall_program = SNPCall_program,
      GenotypeCall_program = GenotypeCall_program,
      CountsFrom = CountsFrom,
      simulated_phases = simulated_phases
  }

  call utilsR.ErrorsReport{
    input:
      onemap_obj = select_onemap_obj,
      simu_onemap_obj = simu_onemap_obj,
      SNPCall_program = SNPCall_program,
      GenotypeCall_program = GenotypeCall_program,
      CountsFrom = CountsFrom
  }

  output{
    File RDatas = MapsReport.maps_RData
    File maps_report = MapsReport.maps_report
    File times = MapsReport.times
    File filters_report = FiltersReport.filters_report
    File errors_report = ErrorsReport.errors_report
  }
}

task GQProbs{
  input{
    File vcf_file
    File onemap_obj
    String cross
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

      onemap_obj <- load("~{onemap_obj}")
      onemap_obj <- get(onemap_obj)

      # MAPS REPORT - GQ
      gq <- extract_depth(vcfR.object=vcf,
                               onemap.object=onemap_obj,
                               vcf.par="GQ",
                               parent1="P1",
                               parent2="P2",
                               f1 = f1,
                               recovering=FALSE)

      gq_onemap_obj <- create_probs(onemap.obj = onemap_obj, genotypes_errors=gq)
      save(gq_onemap_obj, file="gq_onemap_obj.RData")

    RSCRIPT

  >>>
  runtime{
    docker:"cristaniguti/onemap_workflows"
    time:"10:00:00"
    mem:"50GB"
    cpu:1
  }

  output{
    File gq_onemap_obj = "gq_onemap_obj.RData"
  }
}
