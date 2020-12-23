version 1.0

import "./utilsR.wdl" as utilsR


workflow SnpBasedGenotypingSimulatedMaps {

  input {
    File simu_onemap_obj
    File onemap_obj
    File vcf_file
    File ref_alt_alleles
    File simulated_phases
    String SNPCall_program
    String genotyping_program
    String CountsFrom
    String cross
    File? multi_obj
  }

  call OnemapProbsSimulated {
    input:
      method=genotyping_program,
      vcf_file = vcf_file,
      onemap_obj = onemap_obj,
      cross = cross
  }

  call utilsR.GlobalError{
    input:
      onemap_obj = OnemapProbsSimulated.onemap_obj_out
  }

  Array[String] methods                         = [genotyping_program, genotyping_program + "0.05"]
  Array[File] objects                           = [OnemapProbsSimulated.onemap_obj_out, GlobalError.error_onemap_obj]
  Array[Pair[String, File]] methods_and_objects = zip(methods, objects)

  scatter (item in methods_and_objects) {
  
       if (defined(multi_obj)) {
           call utilsR.AddMultiallelics{
             input:
               onemap_obj_multi = multi_obj,
               onemap_obj_bi = item.right
           }
       }
        
       File select_onemap_obj = select_first([AddMultiallelics.onemap_obj_both, item.right])  

       call utilsR.FiltersReport {
            input:
              onemap_obj = select_onemap_obj,
              SNPCall_program = SNPCall_program,
              GenotypeCall_program = item.left,
              CountsFrom = CountsFrom
        }

        call utilsR.MapsReport {
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
              CountsFrom = CountsFrom
          }

   }

   output {
      Array[File] RDatas = MapsReport.maps_RData
      Array[File] maps_report = MapsReport.maps_report
      Array[File] times = MapsReport.times
      Array[File] filters_report = FiltersReport.filters_report
      Array[File] errors_report = ErrorsReport.errors_report
   }
}


task OnemapProbsSimulated {
  input {
    File vcf_file
    File onemap_obj
    String method
    String cross
  }

  command <<<
     R --vanilla --no-save <<RSCRIPT
       library(onemap)
       library(vcfR)
       library(genotyping4onemap)
       method <- "~{method}"
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

       if (method == "updog") {
            out_onemap_obj <- updog_genotype(vcfR.object=vcf,
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
       } else if (method == "supermassa") {
            out_onemap_obj <- supermassa_genotype(vcfR.object=vcf,
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

       } else if (method == "polyrad") {
            out_onemap_obj <- polyRAD_genotype(vcf="~{vcf_file}",
                                                onemap.obj=onemap_obj,
                                                parent1="P1",
                                                parent2="P2",
                                                f1 = f1,
                                                crosstype= cross,
                                                global_error = NULL,
                                                use_genotypes_errors = FALSE,
                                                use_genotypes_probs = TRUE)

       }

       save(out_onemap_obj, file="~{method}_onemap_obj.RData")
     RSCRIPT

  >>>

  runtime{
    docker:"cristaniguti/onemap_workflows"
    time:"96:00:00"
    mem:"--nodes=1"
    cpu:20
  }

  output {
    File onemap_obj_out = "~{method}_onemap_obj.RData"
  }
}
