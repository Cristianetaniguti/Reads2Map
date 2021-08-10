version 1.0

import "structs/maps_empS.wdl"

import "tasks/default_maps_emp.wdl" as default
import "tasks/snpcaller_maps_emp.wdl" as snpcaller
import "tasks/gusmap_maps_emp.wdl" as gusmap
import "tasks/utils.wdl" as utils
import "tasks/utilsR.wdl" as utilsR

import "tasks/genotyping.wdl" as genotyping

struct PopulationAnalysis {
    String method
    File vcf
    File bam
    File? multi
}

workflow Maps {

    input {
        Dataset dataset
        File gatk_vcf
        File freebayes_vcf
        File gatk_vcf_bam_counts
        File freebayes_vcf_bam_counts
        String? filters
        File? gatk_multi
        File? freebayes_multi
        Int max_cores
    }

    if (defined(filters)) {
        call utils.ApplyRandomFilters {
            input:
                gatk_vcf = gatk_vcf,
                freebayes_vcf = freebayes_vcf,
                gatk_vcf_bam_counts = gatk_vcf_bam_counts,
                freebayes_vcf_bam_counts = freebayes_vcf_bam_counts,
                filters = filters,
                chromosome = dataset.chromosome
        }
    }

    File filtered_gatk_vcf = select_first([ApplyRandomFilters.gatk_vcf_filt, gatk_vcf])
    File filtered_gatk_vcf_bamcounts = select_first([ApplyRandomFilters.gatk_vcf_bam_counts_filt, gatk_vcf_bam_counts])
    File filtered_freebayes_vcf = select_first([ApplyRandomFilters.freebayes_vcf_filt, freebayes_vcf])
    File filtered_freebayes_vcf_bamcounts = select_first([ApplyRandomFilters.freebayes_vcf_bam_counts_filt, freebayes_vcf_bam_counts])

    PopulationAnalysis gatk_processing = {"multi": gatk_multi,"method": "gatk", "vcf": filtered_gatk_vcf, "bam": filtered_gatk_vcf_bamcounts}
    PopulationAnalysis freebayes_processing = {"multi": freebayes_multi,"method": "freebayes", "vcf": filtered_freebayes_vcf, "bam": filtered_freebayes_vcf_bamcounts}

    scatter (analysis in [gatk_processing, freebayes_processing]) {

        call utilsR.vcf2onemap {
            input:
                vcf_file = analysis.vcf,
                cross = dataset.cross,
                SNPCall_program = analysis.method,
                parent1 = dataset.parent1,
                parent2 = dataset.parent2
        }

        call utilsR.MultiVcf2onemap{
             input:
                multi = analysis.multi,
                cross = dataset.cross,
                SNPCall_program = analysis.method,
                parent1 = dataset.parent1,
                parent2 = dataset.parent2,
                multiallelics = dataset.multiallelics
        }

        call default.DefaultMaps {
            input:
                onemap_obj = vcf2onemap.onemap_obj,
                vcfR_obj = vcf2onemap.vcfR_obj,
                parent1 = dataset.parent1,
                parent2 = dataset.parent2,
                SNPCall_program = analysis.method,
                CountsFrom = "vcf",
                chromosome = dataset.chromosome,
                multi_obj = MultiVcf2onemap.onemap_obj,
                multiallelics = dataset.multiallelics,
                max_cores = max_cores
        }

        call snpcaller.SNPCallerMaps {
            input:
                onemap_obj = vcf2onemap.onemap_obj,
                vcf_file = analysis.vcf,
                cross = dataset.cross,
                SNPCall_program = analysis.method,
                GenotypeCall_program = "SNPCaller",
                CountsFrom = "vcf",
                parent1 = dataset.parent1,
                parent2 = dataset.parent2,
                chromosome = dataset.chromosome,
                multi_obj = MultiVcf2onemap.onemap_obj,
                multiallelics = dataset.multiallelics,
                max_cores = max_cores
        }

        Map[String, File] vcfs = {"vcf": analysis.vcf, "bam": analysis.bam}

        scatter (origin in ["vcf", "bam"]) {
            call genotyping.SnpBasedGenotypingMaps as UpdogMaps {
                input:
                    onemap_obj = vcf2onemap.onemap_obj,
                    vcf_file = vcfs[origin],
                    SNPCall_program = analysis.method,
                    GenotypeCall_program = "updog",
                    CountsFrom = origin,
                    cross = dataset.cross,
                    parent1 = dataset.parent1,
                    parent2 = dataset.parent2,
                    chromosome = dataset.chromosome,
                    multi_obj = MultiVcf2onemap.onemap_obj,
                    multiallelics = dataset.multiallelics,
                    max_cores = max_cores
            }

            call genotyping.SnpBasedGenotypingMaps as SupermassaMaps {
                input:
                    onemap_obj = vcf2onemap.onemap_obj,
                    vcf_file = vcfs[origin],
                    SNPCall_program = analysis.method,
                    GenotypeCall_program = "supermassa",
                    CountsFrom = origin,
                    cross = dataset.cross,
                    parent1 = dataset.parent1,
                    parent2 = dataset.parent2,
                    chromosome = dataset.chromosome,
                    multi_obj = MultiVcf2onemap.onemap_obj,
                    multiallelics = dataset.multiallelics,
                    max_cores = max_cores
            }

            call genotyping.SnpBasedGenotypingMaps as PolyradMaps {
                input:
                    onemap_obj = vcf2onemap.onemap_obj,
                    vcf_file = vcfs[origin],
                    SNPCall_program = analysis.method,
                    GenotypeCall_program = "polyrad",
                    CountsFrom = origin,
                    cross = dataset.cross,
                    parent1 = dataset.parent1,
                    parent2 = dataset.parent2,
                    chromosome = dataset.chromosome,
                    multi_obj = MultiVcf2onemap.onemap_obj,
                    multiallelics = dataset.multiallelics,
                    max_cores = max_cores
            }
        }

        call gusmap.GusmapMaps {
            input:
                vcf_file = analysis.vcf,
                new_vcf_file = analysis.bam,
                SNPCall_program = analysis.method,
                GenotypeCall_program = "gusmap",
                parent1 = dataset.parent1,
                parent2 = dataset.parent2,
                max_cores = max_cores
        }
    }

    call JointReports {
        input:
            default_RDatas = flatten(DefaultMaps.RDatas),
            default_maps_report = flatten(DefaultMaps.maps_report),
            default_filters_report = flatten(DefaultMaps.filters_report),
            default_errors_report = flatten(DefaultMaps.errors_report),
            default_times_report = flatten(DefaultMaps.times),
            SNPCaller_RDatas = SNPCallerMaps.RDatas,
            SNPCaller_maps_report = SNPCallerMaps.maps_report,
            SNPCaller_filters_report = SNPCallerMaps.filters_report,
            SNPCaller_errors_report = SNPCallerMaps.errors_report,
            SNPCaller_times_report = SNPCallerMaps.times,
            Updog_RDatas = flatten(flatten(UpdogMaps.RDatas)),
            Updog_maps_report = flatten(flatten(UpdogMaps.maps_report)),
            Updog_filters_report = flatten(flatten(UpdogMaps.filters_report)),
            Updog_errors_report = flatten(flatten(UpdogMaps.errors_report)),
            Updog_times_report = flatten(flatten(UpdogMaps.times)),
            Polyrad_RDatas = flatten(flatten(PolyradMaps.RDatas)),
            Polyrad_maps_report = flatten(flatten(PolyradMaps.maps_report)),
            Polyrad_filters_report = flatten(flatten(PolyradMaps.filters_report)),
            Polyrad_errors_report = flatten(flatten(PolyradMaps.errors_report)),
            Polyrad_times_report = flatten(flatten(PolyradMaps.times)),
            Supermassa_RDatas = flatten(flatten(SupermassaMaps.RDatas)),
            Supermassa_maps_report = flatten(flatten(SupermassaMaps.maps_report)),
            Supermassa_filters_report = flatten(flatten(SupermassaMaps.filters_report)),
            Supermassa_errors_report = flatten(flatten(SupermassaMaps.errors_report)),
            Supermassa_times_report = flatten(flatten(SupermassaMaps.times)),
            Gusmap_RDatas = flatten(GusmapMaps.RDatas),
            Gusmap_maps_report = flatten(GusmapMaps.maps_report),
            Gusmap_times_report = flatten(GusmapMaps.times),
            max_cores = max_cores
    }

    output{
        File EmpiricalReads_results = JointReports.EmpiricalReads_results
    }
}

task JointReports{
  input{
    Array[File] default_RDatas
    Array[File] default_maps_report
    Array[File] default_filters_report
    Array[File] default_errors_report
    Array[File] default_times_report
    Array[File] SNPCaller_RDatas
    Array[File] SNPCaller_maps_report
    Array[File] SNPCaller_filters_report
    Array[File] SNPCaller_errors_report
    Array[File] SNPCaller_times_report
    Array[File] Updog_RDatas
    Array[File] Updog_maps_report
    Array[File] Updog_filters_report
    Array[File] Updog_errors_report
    Array[File] Updog_times_report
    Array[File] Polyrad_RDatas
    Array[File] Polyrad_maps_report
    Array[File] Polyrad_filters_report
    Array[File] Polyrad_errors_report
    Array[File] Polyrad_times_report
    Array[File] Supermassa_RDatas
    Array[File] Supermassa_maps_report
    Array[File] Supermassa_filters_report
    Array[File] Supermassa_errors_report
    Array[File] Supermassa_times_report
    Array[File] Gusmap_RDatas
    Array[File] Gusmap_maps_report
    Array[File] Gusmap_times_report
    Int max_cores
  }

  command <<<
     R --vanilla --no-save <<RSCRIPT

     # packages
      library(tidyr)
      library(stringr)
      library(vroom)
      library(largeList)

      # function
      joint_reports <- function(default, 
                                snpcaller, 
                                updog, 
                                polyrad, 
                                supermassa,
                                gusmap = NULL){

        default    <- str_split(default, ";", simplify = T)
        SNPCaller  <- str_split(snpcaller, ";", simplify = T)
        Updog      <- str_split(updog, ";", simplify = T)
        Polyrad    <- str_split(polyrad, ";", simplify = T)
        Supermassa <- str_split(supermassa, ";", simplify = T)

        if(is.null(gusmap)){
          files <- c(default, SNPCaller, Updog, Polyrad, Supermassa)
        } else {
          Gusmap <- str_split(gusmap, ";", simplify = T)
          files <- c(default, SNPCaller, Updog, Polyrad, Supermassa, Gusmap)
        }

        joint <- vroom(files, num_threads = ~{max_cores})
        return(joint)
      }
    
      maps_report <- joint_reports(default = "~{sep=";" default_maps_report}", 
                                  snpcaller = "~{sep=";" SNPCaller_maps_report}", 
                                  updog = "~{sep=";" Updog_maps_report}", 
                                  polyrad = "~{sep=";" Polyrad_maps_report}", 
                                  supermassa = "~{sep=";" Supermassa_maps_report}",
                                  gusmap = "~{sep=";" Gusmap_maps_report}")

      filters_report <- joint_reports(default = "~{sep=";" default_filters_report}", 
                            snpcaller = "~{sep=";" SNPCaller_filters_report}", 
                            updog = "~{sep=";" Updog_filters_report}", 
                            polyrad = "~{sep=";" Polyrad_filters_report}", 
                            supermassa = "~{sep=";" Supermassa_filters_report}")

      errors_report <- joint_reports(default = "~{sep=";" default_errors_report}", 
                                     snpcaller = "~{sep=";" SNPCaller_errors_report}", 
                                     updog = "~{sep=";" Updog_errors_report}", 
                                     polyrad = "~{sep=";" Polyrad_errors_report}", 
                                     supermassa = "~{sep=";" Supermassa_errors_report}")

 
       times_report <- joint_reports(default = "~{sep=";" default_times_report}", 
                                   snpcaller = "~{sep=";" SNPCaller_times_report}", 
                                   updog = "~{sep=";" Updog_times_report}", 
                                   polyrad = "~{sep=";" Polyrad_times_report}", 
                                   supermassa = "~{sep=";" Supermassa_times_report}",
                                   gusmap = "~{sep=";" Gusmap_times_report}")

     # RDatas need to be load
      default    <- str_split("~{sep=";" default_RDatas}", ";", simplify = T)
      SNPCaller  <- str_split("~{sep=";" SNPCaller_RDatas}", ";", simplify = T)
      Updog      <- str_split("~{sep=";" Updog_RDatas}", ";", simplify = T)
      Polyrad    <- str_split("~{sep=";" Polyrad_RDatas}", ";", simplify = T)
      Supermassa <- str_split("~{sep=";" Supermassa_RDatas}", ";", simplify = T)
      Gusmap <- str_split("~{sep=";" Gusmap_RDatas}", ";", simplify = T)

      RDatas_names <- c(default, SNPCaller, Updog, Polyrad, Supermassa, Gusmap)

      all_RDatas <- list()
      for(i in 1:length(RDatas_names)){
         map_temp <- load(RDatas_names[i])
         all_RDatas[[i]] <- get(map_temp)
      }

      names(all_RDatas) <- sapply(RDatas_names, basename)
      gusmap_RDatas <- all_RDatas[grep("gusmap", names(all_RDatas))]
      RDatas <- all_RDatas[-grep("gusmap", names(all_RDatas))]

      # Converting OneMap sequencig objects to list. LargeList do not accept other class
      # Also because of this gusmap is separated, because the developers worked with enviroments, not classes

      for(i in 1:length(RDatas)){
        class(RDatas[[i]]) <- "list"
      }

      saveList(RDatas, file = "sequences_emp.llo", append=FALSE, compress=TRUE)

      new_names <- names(all_RDatas)
      vroom_write(as.data.frame(new_names), "names.tsv.gz")
      save(gusmap_RDatas, file = "gusmap_RDatas.RData")

      # Outputs
      vroom_write(errors_report, "data1_depths_geno_prob.tsv.gz", num_threads = ~{max_cores})
      vroom_write(maps_report, "data2_maps.tsv.gz", num_threads = ~{max_cores})
      vroom_write(filters_report, "data3_filters.tsv.gz", num_threads = ~{max_cores})
      vroom_write(times_report, "data4_times.tsv.gz", num_threads = ~{max_cores})
      
      system("mkdir EmpiricalReads_results")
      system("mv gusmap_RDatas.RData sequences_emp.llo data1_depths_geno_prob.tsv.gz data2_maps.tsv.gz data3_filters.tsv.gz data4_times.tsv.gz names.tsv.gz EmpiricalReads_results")
      system("tar -czvf EmpiricalReads_results.tar.gz EmpiricalReads_results")

     RSCRIPT

  >>>

  runtime{
    docker:"cristaniguti/reads2map:0.0.1"
    time:"10:00:00"
    mem:"80GB"
    cpu:1
  }

  output{
    File EmpiricalReads_results = "EmpiricalReads_results.tar.gz"
  }
}
