version 1.0

import "../structs/maps_empS.wdl"

import "./snpcaller_maps_emp.wdl" as snpcaller
import "./gusmap_maps_emp.wdl" as gusmap
import "./utils.wdl" as utils
import "./utilsR.wdl" as utilsR

import "./genotyping.wdl" as genotyping

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
        Int max_cores
        File merged_bam
        File reference
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

    PopulationAnalysis gatk_processing = {"method": "gatk", "vcf": filtered_gatk_vcf, "bam": filtered_gatk_vcf_bamcounts}
    PopulationAnalysis freebayes_processing = {"method": "freebayes", "vcf": filtered_freebayes_vcf, "bam": filtered_freebayes_vcf_bamcounts}

    # Re-Genotyping with updog, supermassa and polyrad; and building maps with onemap
    scatter (analysis in [gatk_processing, freebayes_processing]) {

        Map[String, File] vcfs = {"vcf": analysis.vcf, "bam": analysis.bam}

        scatter (origin in ["vcf", "bam"]) {

            call genotyping.onemapMaps as updogMaps {
                input:
                    vcf_file = vcfs[origin],
                    SNPCall_program = analysis.method,
                    GenotypeCall_program = "updog",
                    CountsFrom = origin,
                    cross = dataset.cross,
                    parent1 = dataset.parent1,
                    parent2 = dataset.parent2,
                    chromosome = dataset.chromosome,
                    multiallelics = dataset.multiallelics,
                    max_cores = max_cores,
                    reference = reference,
                    merged_bam = merged_bam
            }

            call genotyping.onemapMaps as supermassaMaps {
                input:
                    vcf_file = vcfs[origin],
                    SNPCall_program = analysis.method,
                    GenotypeCall_program = "supermassa",
                    CountsFrom = origin,
                    cross = dataset.cross,
                    parent1 = dataset.parent1,
                    parent2 = dataset.parent2,
                    chromosome = dataset.chromosome,
                    multiallelics = dataset.multiallelics,
                    max_cores = max_cores,
                    reference = reference,
                    merged_bam = merged_bam
            }

            call genotyping.onemapMaps as polyradMaps {
                input:
                    vcf_file = vcfs[origin],
                    SNPCall_program = analysis.method,
                    GenotypeCall_program = "polyrad",
                    CountsFrom = origin,
                    cross = dataset.cross,
                    parent1 = dataset.parent1,
                    parent2 = dataset.parent2,
                    chromosome = dataset.chromosome,
                    multiallelics = dataset.multiallelics,
                    max_cores = max_cores,
                    reference = reference,
                    merged_bam = merged_bam
            }
        }

       # Build maps with GUSMap
       call gusmap.gusmapMaps {
            input:
              vcf_file = analysis.vcf,
              new_vcf_file = analysis.bam,
              SNPCall_program = analysis.method,
              GenotypeCall_program = "gusmap",
              parent1 = dataset.parent1,
              parent2 = dataset.parent2,
              max_cores = max_cores
        }

        call snpcaller.SNPCallerMaps {
            input:
                vcf_file = analysis.vcf,
                cross = dataset.cross,
                SNPCall_program = analysis.method,
                GenotypeCall_program = "SNPCaller",
                CountsFrom = "vcf",
                parent1 = dataset.parent1,
                parent2 = dataset.parent2,
                chromosome = dataset.chromosome,
                multiallelics = dataset.multiallelics,
                max_cores = max_cores,
                reference = reference,
                merged_bam = merged_bam
        }
    }

    # Compress files
    call JointReports {
        input:
            SNPCaller = SNPCallerMaps.tar_gz_report,
            updog = flatten(updogMaps.tar_gz_report),
            polyrad = flatten(polyradMaps.tar_gz_report),
            supermassa = flatten(supermassaMaps.tar_gz_report),
            gusmap = gusmapMaps.tar_gz_report,
            max_cores = max_cores
    }

    output{
        File EmpiricalReads_results = JointReports.EmpiricalReads_results
    }
}

task JointReports{
  input{
    Array[File] SNPCaller
    Array[File] updog
    Array[File] polyrad
    Array[File] supermassa
    Array[File] gusmap
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
        updog      <- str_split(updog, ";", simplify = T)
        polyrad    <- str_split(polyrad, ";", simplify = T)
        supermassa <- str_split(supermassa, ";", simplify = T)

        if(is.null(gusmap)){
          files <- c(default, SNPCaller, updog, polyrad, supermassa)
        } else {
          gusmap <- str_split(gusmap, ";", simplify = T)
          files <- c(default, SNPCaller, updog, polyrad, supermassa, gusmap)
        }

        joint <- vroom(files, num_threads = ~{max_cores})
        return(joint)
      }
    
      maps_report <- joint_reports(snpcaller = "~{sep=";" SNPCaller}", 
                                  updog = "~{sep=";" updog}", 
                                  polyrad = "~{sep=";" polyrad}", 
                                  supermassa = "~{sep=";" supermassa}",
                                  gusmap = "~{sep=";" gusmap}")

    #   filters_report <- joint_reports(snpcaller = "~{sep=";" SNPCaller}", 
    #                         updog = "~{sep=";" updog}", 
    #                         polyrad = "~{sep=";" polyrad}", 
    #                         supermassa = "~{sep=";" supermassa}")

    #   errors_report <- joint_reports(snpcaller = "~{sep=";" SNPCaller}", 
    #                                  updog = "~{sep=";" updog}", 
    #                                  polyrad = "~{sep=";" polyrad}", 
    #                                  supermassa = "~{sep=";" supermassa}")

 
    #    times_report <- joint_reports(snpcaller = "~{sep=";" SNPCaller}", 
    #                                updog = "~{sep=";" updog}", 
    #                                polyrad = "~{sep=";" polyrad}", 
    #                                supermassa = "~{sep=";" supermassa}",
    #                                gusmap = "~{sep=";" gusmap}")

    #  # RDatas need to be load
    #   SNPCaller  <- str_split("~{sep=";" SNPCaller}", ";", simplify = T)
    #   updog      <- str_split("~{sep=";" updog}", ";", simplify = T)
    #   polyrad    <- str_split("~{sep=";" polyrad}", ";", simplify = T)
    #   supermassa <- str_split("~{sep=";" supermassa}", ";", simplify = T)
    #   gusmap <- str_split("~{sep=";" gusmap}", ";", simplify = T)

    #   RDatas_names <- c(default, SNPCaller, updog, polyrad, supermassa, gusmap)

    #   all <- list()
    #   for(i in 1:length(RDatas_names)){
    #      map_temp <- load(RDatas_names[i])
    #      all_RDatas[[i]] <- get(map_temp)
    #   }

    #   names(all_RDatas) <- sapply(RDatas_names, basename)
    #   gusmap_RDatas <- all_RDatas[grep("gusmap", names(all_RDatas))]
    #   RDatas <- all_RDatas[-grep("gusmap", names(all_RDatas))]

    #   # Converting onemap sequencig objects to list. LargeList do not accept other class
    #   # Also because of this gusmap is separated, because the developers worked with enviroments, not classes

    #   for(i in 1:length(RDatas)){
    #     class(RDatas[[i]]) <- "list"
    #   }

    #   saveList(RDatas, file = "sequences_emp.llo", append=FALSE, compress=TRUE)

    #   new_names <- names(all_RDatas)
    #   vroom_write(as.data.frame(new_names), "names.tsv.gz")
    #   save(gusmap_RDatas, file = "gusmap_RDatas.RData")

    #   # Outputs
    #   vroom_write(errors_report, "data1_depths_geno_prob.tsv.gz", num_threads = ~{max_cores})
    #   vroom_write(maps_report, "data2_maps.tsv.gz", num_threads = ~{max_cores})
    #   vroom_write(filters_report, "data3_filters.tsv.gz", num_threads = ~{max_cores})
    #   vroom_write(times_report, "data4_times.tsv.gz", num_threads = ~{max_cores})
      
    #   system("mkdir EmpiricalReads_results")
    #   system("mv gusmap_RDatas.RData sequences_emp.llo data1_depths_geno_prob.tsv.gz data2_maps.tsv.gz data3_filters.tsv.gz data4_times.tsv.gz names.tsv.gz EmpiricalReads_results")
    #   system("tar -czvf EmpiricalReads_results.tar.gz EmpiricalReads_results")

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
