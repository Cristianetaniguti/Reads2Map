version 1.0

import "../structs/maps_empS.wdl"

import "default_maps_emp.wdl" as default
import "snpcaller_maps_emp.wdl" as snpcaller
import "updog_maps_emp.wdl" as updog
import "polyrad_maps_emp.wdl" as polyrad
import "supermassa_maps_emp.wdl" as supermassa
import "gusmap_maps_emp.wdl" as gusmap
import "utils.wdl" as utils
import "utilsR.wdl" as utilsR


struct PopulationAnalysis {
    String method
    File vcf
    File bam
}

workflow Maps {

    input {
        Dataset dataset
        File gatk_vcf
        File freebayes_vcf
        File gatk_vcf_bam_counts
        File freebayes_vcf_bam_counts
        String? filters
    }

    if (defined(filters)) {
        call utils.ApplyRandomFilters {
            input:
                gatk_vcf = gatk_vcf,
                freebayes_vcf = freebayes_vcf,
                gatk_vcf_bam_counts = gatk_vcf_bam_counts,
                freebayes_vcf_bam_counts = freebayes_vcf_bam_counts,
                filters = filters
        }
    }

    File filtered_gatk_vcf = select_first([ApplyRandomFilters.gatk_vcf_filt, gatk_vcf])
    File filtered_gatk_vcf_bamcounts = select_first([ApplyRandomFilters.gatk_vcf_bam_counts_filt, gatk_vcf_bam_counts])
    File filtered_freebayes_vcf = select_first([ApplyRandomFilters.freebayes_vcf_filt, freebayes_vcf])
    File filtered_freebayes_vcf_bamcounts = select_first([ApplyRandomFilters.freebayes_vcf_bam_counts_filt, freebayes_vcf_bam_counts])

    PopulationAnalysis gatk_processing = {"method": "gatk", "vcf": filtered_gatk_vcf, "bam": filtered_gatk_vcf_bamcounts}
    PopulationAnalysis freebayes_processing = {"method": "freebayes", "vcf": filtered_freebayes_vcf, "bam": filtered_freebayes_vcf_bamcounts}

    scatter (analysis in [gatk_processing, freebayes_processing]) {

        call utilsR.vcf2onemap {
            input:
                vcf_file = analysis.vcf,
                cross = dataset.cross,
                SNPCall_program = analysis.method,
                parent1 = dataset.parent1,
                parent2 = dataset.parent2
        }

        call default.DefaultMaps {
            input:
                onemap_obj = vcf2onemap.onemap_obj,
                vcfR_obj = vcf2onemap.vcfR_obj,
                parent1 = dataset.parent1,
                parent2 = dataset.parent2,
                SNPCall_program = analysis.method,
                GenotypeCall_program = "default",
                CountsFrom = "vcf",
                cross = dataset.cross,
                chromosome = dataset.chromosome
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
                chromosome = dataset.chromosome
        }

        Map[String, File] vcfs = {"vcf": analysis.vcf, "bam": analysis.bam}

        scatter (origin in ["vcf", "bam"]) {
            call updog.UpdogMaps {
                input:
                    onemap_obj = vcf2onemap.onemap_obj,
                    vcf_file = vcfs[origin],
                    SNPCall_program = analysis.method,
                    GenotypeCall_program = "updog",
                    CountsFrom = origin,
                    cross = dataset.cross,
                    parent1 = dataset.parent1,
                    parent2 = dataset.parent2,
                    chromosome = dataset.chromosome
            }

            call supermassa.SupermassaMaps {
                input:
                    onemap_obj = vcf2onemap.onemap_obj,
                    vcf_file = vcfs[origin],
                    SNPCall_program = analysis.method,
                    GenotypeCall_program = "supermassa",
                    CountsFrom = origin,
                    cross = dataset.cross,
                    parent1 = dataset.parent1,
                    parent2 = dataset.parent2,
                    chromosome = dataset.chromosome
            }

            call polyrad.PolyradMaps {
                input:
                    onemap_obj = vcf2onemap.onemap_obj,
                    vcf_file = vcfs[origin],
                    SNPCall_program = analysis.method,
                    GenotypeCall_program = "polyrad",
                    CountsFrom = origin,
                    cross = dataset.cross,
                    parent1 = dataset.parent1,
                    parent2 = dataset.parent2,
                    chromosome = dataset.chromosome
            }
        }

        call gusmap.GusmapMaps {
            input:
                vcf_file = analysis.vcf,
                new_vcf_file = analysis.bam,
                SNPCall_program = analysis.method,
                GenotypeCall_program = "gusmap",
                parent1 = dataset.parent1,
                parent2 = dataset.parent2
        }
    }

    call JointReports {
        input:
            default_RDatas = flatten(DefaultMaps.RDatas),
            default_maps_report = flatten(DefaultMaps.maps_report),
            default_filters_report = flatten(DefaultMaps.filters_report),
            default_errors_report = flatten(DefaultMaps.errors_report),
            default_times = flatten(DefaultMaps.times),
            SNPCaller_RDatas = SNPCallerMaps.RDatas,
            SNPCaller_maps_report = SNPCallerMaps.maps_report,
            SNPCaller_filters_report = SNPCallerMaps.filters_report,
            SNPCaller_errors_report = SNPCallerMaps.errors_report,
            SNPCaller_times = SNPCallerMaps.times,
            Updog_RDatas = flatten(flatten(UpdogMaps.RDatas)),
            Updog_maps_report = flatten(flatten(UpdogMaps.maps_report)),
            Updog_filters_report = flatten(flatten(UpdogMaps.filters_report)),
            Updog_errors_report = flatten(flatten(UpdogMaps.errors_report)),
            Updog_times = flatten(flatten(UpdogMaps.times)),
            Polyrad_RDatas = flatten(flatten(PolyradMaps.RDatas)),
            Polyrad_maps_report = flatten(flatten(PolyradMaps.maps_report)),
            Polyrad_filters_report = flatten(flatten(PolyradMaps.filters_report)),
            Polyrad_errors_report = flatten(flatten(PolyradMaps.errors_report)),
            Polyrad_times = flatten(flatten(PolyradMaps.times)),
            Supermassa_RDatas = flatten(flatten(SupermassaMaps.RDatas)),
            Supermassa_maps_report = flatten(flatten(SupermassaMaps.maps_report)),
            Supermassa_filters_report = flatten(flatten(SupermassaMaps.filters_report)),
            Supermassa_errors_report = flatten(flatten(SupermassaMaps.errors_report)),
            Supermassa_times = flatten(flatten(SupermassaMaps.times)),
            Gusmap_RDatas = flatten(GusmapMaps.RDatas),
            Gusmap_maps_report = flatten(GusmapMaps.maps_report),
            Gusmap_times = flatten(GusmapMaps.times)
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
    Array[File] default_times
    Array[File] SNPCaller_RDatas
    Array[File] SNPCaller_maps_report
    Array[File] SNPCaller_filters_report
    Array[File] SNPCaller_errors_report
    Array[File] SNPCaller_times
    Array[File] Updog_RDatas
    Array[File] Updog_maps_report
    Array[File] Updog_filters_report
    Array[File] Updog_errors_report
    Array[File] Updog_times
    Array[File] Polyrad_RDatas
    Array[File] Polyrad_maps_report
    Array[File] Polyrad_filters_report
    Array[File] Polyrad_errors_report
    Array[File] Polyrad_times
    Array[File] Supermassa_RDatas
    Array[File] Supermassa_maps_report
    Array[File] Supermassa_filters_report
    Array[File] Supermassa_errors_report
    Array[File] Supermassa_times
    Array[File] Gusmap_RDatas
    Array[File] Gusmap_maps_report
    Array[File] Gusmap_times
  }

  command <<<
     R --vanilla --no-save <<RSCRIPT
      # I needed to split in groups because of R character limit size
      system("cat ~{sep= ' ' default_maps_report} ~{sep= ' ' SNPCaller_maps_report} > temp_map1")
      system("cat ~{sep= ' ' Updog_maps_report} ~{sep= ' ' Polyrad_maps_report} > temp_map2")
      system("cat ~{sep= ' ' Supermassa_maps_report}  ~{sep= ' ' Gusmap_maps_report} > temp_map3")
      system("cat temp_map1 temp_map2  temp_map3 > all_maps.txt")

      system("cat ~{sep= ' ' default_filters_report} ~{sep= ' ' SNPCaller_filters_report} > temp_filters1")
      system("cat ~{sep= ' ' Updog_filters_report} ~{sep= ' ' Polyrad_filters_report} > temp_filters2")
      system("cat ~{sep= ' ' Supermassa_filters_report} > temp_filters3")
      system("cat temp_filters1 temp_filters2 temp_filters3 > all_filters.txt")

      system("cat ~{sep= ' '  default_errors_report} ~{sep= ' ' SNPCaller_errors_report} > temp_errors1")
      system("cat ~{sep= ' ' Updog_errors_report} ~{sep= ' ' Polyrad_errors_report} > temp_errors2")
      system("cat ~{sep= ' ' Supermassa_errors_report} > temp_errors3")
      system("cat temp_errors1 temp_errors2 temp_errors3 > all_errors.txt")

      system("cat ~{sep= ' '  default_times} ~{sep= ' ' SNPCaller_times} > temp_times1")
      system("cat ~{sep= ' ' Updog_times} ~{sep= ' ' Polyrad_times} > temp_times2")
      system("cat ~{sep= ' ' Supermassa_times}  ~{sep= ' ' Gusmap_times} > temp_times3")
      system("cat temp_times1 temp_times2 temp_times3 > all_times.txt")

      # RDatas need to be load
      system("cp ~{sep= ' ' default_RDatas} ~{sep= ' ' SNPCaller_RDatas}  .")
      system("cp ~{sep= ' ' Updog_RDatas}  ~{sep= ' ' Polyrad_RDatas} .")
      system("cp ~{sep= ' ' Supermassa_RDatas} ~{sep= ' ' Gusmap_RDatas} .")

     library(tidyr)
     library(largeList)
     library(data.table)
     source("/opt/scripts/functions_empirical.R")

     Genocall <- c("default", "SNPCaller", "updog", "supermassa", "polyrad", "gusmap",
                   "default0.05", "updog0.05", "supermassa0.05", "polyrad0.05")
     fake <- c(TRUE, FALSE)
     CountsFrom <- c("vcf", "bam")
     SNPCall <- c("gatk", "freebayes")

      df <- data.frame(SNPCall, CountsFrom, Genocall)
      df <- tidyr::expand(df, SNPCall, CountsFrom, Genocall)
      df <- as.data.frame(df)
      df <- df[-which((df[,3] == "default" | df[,3] == "default0.05" | df[,3] == "SNPCaller" ) & df[,2] == "bam"),]
      RDatas_names <- paste0("map_",df[,1],"_",df[,2], "_", df[,3],".RData")

      all_RDatas <- list()
      for(i in 1:length(RDatas_names)){
         map_temp <- load(RDatas_names[i])
         all_RDatas[[i]] <- get(map_temp)
      }

      names(all_RDatas) <- RDatas_names
      gusmap_RDatas <- all_RDatas[grep("gusmap", names(all_RDatas))]
      RDatas <- all_RDatas[-grep("gusmap", names(all_RDatas))]

      # Converting OneMap sequencig objects to list. LargeList do not accept other class
      # Also because of this gusmap is separated, because the developers worked with enviroments, not classes

      for(i in 1:length(RDatas)){
        class(RDatas[[i]]) <- "list"
      }

       saveList(RDatas, file = "sequences_emp.llo", append=FALSE, compress=TRUE)

       new_names <- names(all_RDatas)
       saveRDS(new_names, file = "names.rds")
       save(gusmap_RDatas, file = "gusmap_RDatas.RData")

       all_errors <- fread("all_errors.txt")
       colnames(all_errors) <- c("SNPCall", "CountsFrom", "GenoCall", "mks", "ind", "alt", "ref",
                                 "gt.onemap", "gt.vcf", "A", "AB", "BA", "B")
       all_errors <- fix_genocall_names(all_errors)
       saveRDS(all_errors, "data1_depths_geno_prob.rds")

       all_filters <- fread("all_filters.txt")
       colnames(all_filters) <- c("CountsFrom", "SNPCall", "GenoCall",
                                  "miss", "n_markers", "n_markers_selected_chr",
                                  "selected_chr_no_dist", "distorted_markers",
                                  "redundant_markers", "non-grouped_markers")
       all_filters <- fix_genocall_names(all_filters)
       saveRDS(all_filters, "data3_filters.rds")

       all_maps <- fread("all_maps.txt")
       colnames(all_maps) <- c("CountsFrom", "SNPCall", "GenoCall", "mks", "pos", "cm", "mk.type", "phase")
       all_maps <- fix_genocall_names(all_maps)
       saveRDS(all_maps, "data2_maps.rds")

       all_times <- fread("all_times.txt")
       colnames(all_times) <- c("SNPCall", "CountsFrom", "GenoCall", "time")
       all_times <- fix_genocall_names(all_times)
       saveRDS(all_times, "data4_times.rds")

       system("mkdir EmpiricalReads_results")
       system("mv gusmap_RDatas.RData sequences_emp.llo data1_depths_geno_prob.rds data2_maps.rds data3_filters.rds data4.rds data4_times.rds names.rds EmpiricalReads_results")
       system("tar -czvf EmpiricalReads_results.tar.gz EmpiricalReads_results")

     RSCRIPT

  >>>

  runtime{
    docker:"gcr.io/taniguti-backups/onemap:v1"
    time:"48:00:00"
    # mem:"--nodes=1"
    cpu:1
  }

  output{
    File EmpiricalReads_results = "EmpiricalReads_results.tar.gz"
  }
}
