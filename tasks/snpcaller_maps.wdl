version 1.0

import "./utilsR.wdl" as utilsR

workflow SNPCallerMaps{
  input {
     File simu_onemap_obj
     File simu_vcfR
     File onemap_obj
     File vcf_file
     File ref_alt_alleles
     File simulated_phases
     String cross
     String SNPCall_program
     String GenotypeCall_program
     String CountsFrom
     File? multi_obj
     Int seed
     Int depth
     Int max_cores
     String multiallelics
    }


  call SNPCallerProbs {
    input:
      vcf_file = vcf_file,
      onemap_obj = onemap_obj,
      cross = cross
  }

  if (multiallelics == "TRUE") {
      call utilsR.AddMultiallelics{
          input:
            onemap_obj_multi = multi_obj,
            onemap_obj_bi = SNPCallerProbs.probs_onemap_obj
      }
  }

  File select_onemap_obj = select_first([AddMultiallelics.onemap_obj_both, SNPCallerProbs.probs_onemap_obj])

  call utilsR.FiltersReport{
    input:
      onemap_obj = select_onemap_obj,
      SNPCall_program = SNPCall_program,
      GenotypeCall_program = GenotypeCall_program,
      CountsFrom = "vcf",
      seed = seed,
      depth = depth
  }

  call utilsR.MapsReport{
    input:
      onemap_obj = FiltersReport.onemap_obj_filtered,
      ref_alt_alleles = ref_alt_alleles,
      simu_onemap_obj = simu_onemap_obj,
      SNPCall_program = SNPCall_program,
      GenotypeCall_program = GenotypeCall_program,
      CountsFrom = CountsFrom,
      simulated_phases = simulated_phases,
      seed = seed,
      depth = depth,
      max_cores = max_cores
  }

  call utilsR.ErrorsReport{
    input:
      onemap_obj = select_onemap_obj,
      simu_onemap_obj = simu_onemap_obj,
      SNPCall_program = SNPCall_program,
      GenotypeCall_program = GenotypeCall_program,
      CountsFrom = CountsFrom,
      simu_vcfR = simu_vcfR,
      vcfR_obj = SNPCallerProbs.vcfR_obj,
      seed = seed,
      depth = depth,
      max_cores = max_cores
  }

  output {
    File RDatas = MapsReport.maps_RData
    File maps_report = MapsReport.maps_report
    File times = MapsReport.times
    File filters_report = FiltersReport.filters_report
    File errors_report = ErrorsReport.errors_report
  }
}

task SNPCallerProbs{
  input {
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
      save(vcf, file="vcfR_obj.RData")

      onemap_obj <- load("~{onemap_obj}")
      onemap_obj <- get(onemap_obj)

      if(any(grepl("freeBayes", vcf@meta))) par <- "GL" else par <- "PL"

      probs <- extract_depth(vcfR.object=vcf,
                               onemap.object=onemap_obj,
                               vcf.par=par,
                               parent1="P1",
                               parent2="P2",
                               f1 = f1,
                               recovering=FALSE)

      probs_onemap_obj <- create_probs(onemap.obj = onemap_obj, genotypes_probs = probs)
      save(probs_onemap_obj, file="probs_onemap_obj.RData")

    RSCRIPT

  >>>
  runtime {
    docker: "cristaniguti/reads2map"
    preemptible: 3
    memory: "3 GB"
    cpu: 1
  }

  output {
    File probs_onemap_obj = "probs_onemap_obj.RData"
    File vcfR_obj = "vcfR_obj.RData"
  }
}
