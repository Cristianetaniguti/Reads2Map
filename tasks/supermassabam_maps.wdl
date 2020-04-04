version 1.0

import "./utilsR.wdl" as utilsR

workflow SupermassaBamMaps{
  input{
    File simu_onemap_obj
    File vcfR_obj
    File onemap_obj
    File tot_mks
    File real_phases
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    String cMbyMb
    File freebayes_ref_bam
    File freebayes_alt_bam
    File gatk_ref_bam
    File gatk_alt_bam
    String cross
  }
  
  call SupermassaBamProbs{
    input:
      vcfR_obj = vcfR_obj,
      onemap_obj = onemap_obj,
      SNPCall_program = SNPCall_program,
      freebayes_ref_bam = freebayes_ref_bam,
      freebayes_alt_bam = freebayes_alt_bam,
      gatk_ref_bam = gatk_ref_bam,
      gatk_alt_bam = gatk_alt_bam,
      cross=cross
  }
  
  call utilsR.GlobalError{
    input:
      onemap_obj = onemap_obj,
      SNPCall_program = SNPCall_program,
      GenotypeCall_program = GenotypeCall_program,
      CountsFrom = CountsFrom
  }

  Array[String] methods                         = ["supermassa", "supermassa0.05"]
  Array[File] objects                           = [onemap_obj, GlobalError.error_onemap_obj]
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

task SupermassaBamProbs{
  input{
    File vcfR_obj
    File onemap_obj
    String SNPCall_program
    File freebayes_ref_bam
    File freebayes_alt_bam
    File gatk_ref_bam
    File gatk_alt_bam
    String cross
  }
  
  command <<<
     R --vanilla --no-save <<RSCRIPT
       library(onemap)
       library(supermassa4onemap)

       system("cp ~{freebayes_ref_bam} .")
       system("cp ~{freebayes_alt_bam} .")
       system("cp ~{gatk_ref_bam} .")
       system("cp ~{gatk_alt_bam} .")
       
       cross <- "~{cross}"
          
       if(cross == "F1"){
          cross <- "outcross"
          f1 = NULL
       } else if (cross == "F2"){
          cross <- "f2 intercross"
          f1 = "F1"
       }
       
       ## Depths from bam
       depths.alt <- read.table(paste0("~{SNPCall_program}", "_alt_depth_bam.txt"), header = T)
       depths.ref <- read.table(paste0("~{SNPCall_program}", "_ref_depth_bam.txt"), header = T)

       depths <- list("ref" = depths.ref, "alt"=depths.alt)
        
       vcf_temp <- load("~{vcfR_obj}")
       vcf <- get(vcf_temp)
       
       onemap_obj_temp <- load("~{onemap_obj}")
       onemap_obj <- get(onemap_obj_temp)
       
       supermassa_bam_onemap_obj <- supermassa_genotype(vcfR.object=vcf,
                                    onemap.object=onemap_obj,
                                    vcf.par="AD",
                                    parent1="P1",
                                    parent2="P2",
                                    f1 = f1,
                                    recovering=TRUE,
                                    mean_phred=20,
                                    cores=3,
                                    depths=depths,
                                    global_error = NULL,
                                    use_genotypes_errors = FALSE,
                                    use_genotypes_probs = TRUE)
       
       save(supermassa_bam_onemap_obj, file="supermassa_bam_onemap_obj.RData")
  
     RSCRIPT
  >>>
  
  runtime{
    docker:"taniguti/onemap"
  }
  
  output{
    File supermassa_bam_onemap_obj = "supermassa_bam_onemap_obj.RData"
  }
}