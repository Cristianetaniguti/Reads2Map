version 1.0


import "./create_alignment_from_read_simulations.wdl" as simulation
import "./gatk_genotyping.wdl" as gatk
import "./freebayes_genotyping.wdl" as freebayes
import "./utils.wdl" as utils
import "./utilsR.wdl" as utilsR
import "./snpcaller_maps-simulated.wdl" as snpcaller
import "./gusmap_maps-simulated.wdl" as gusmap
import "./genotyping-simulated.wdl" as genotyping


struct PopulationAnalysis {
    String method
    File vcf
    File bam
}

workflow SimulatedSingleFamily {

  input {
    Reference references
    Family family
    Sequencing sequencing
    String? filters
    Int max_cores
    Int chunk_size
    Int ploidy
  }

  call simulation.CreateAlignmentFromSimulation {
    input:
      sequencing = sequencing,
      references=references,
      family=family,
      max_cores = max_cores
  }

  call gatk.GatkGenotyping {
    input:
      bams=CreateAlignmentFromSimulation.bam,
      bais=CreateAlignmentFromSimulation.bai,
      references=references,
      program="gatk",
      vcf_simu = CreateAlignmentFromSimulation.true_vcf,
      seed    = family.seed,
      depth   = sequencing.depth,
      chunk_size = chunk_size,
      ploidy = ploidy
  }

  call freebayes.FreebayesGenotyping {
    input:
      bams=CreateAlignmentFromSimulation.bam,
      bais=CreateAlignmentFromSimulation.bai,
      references=references,
      program="freebayes",
      max_cores = max_cores,
      vcf_simu = CreateAlignmentFromSimulation.true_vcf
  }

  call utilsR.vcf2onemap as truth_vcf {
    input:
      vcf_file = CreateAlignmentFromSimulation.true_vcf,
      cross = family.cross,
      SNPCall_program = "simu",
      parent1 = "P1",
      parent2 = "P2"
  }

  if (defined(filters)) {
      call utils.ApplyRandomFilters {
          input:
              gatk_vcf = GatkGenotyping.vcf_norm,
              freebayes_vcf = FreebayesGenotyping.vcf_norm,
              gatk_vcf_bam_counts = GatkGenotyping.vcf_norm_bamcounts,
              freebayes_vcf_bam_counts = FreebayesGenotyping.vcf_norm_bamcounts,
              filters = filters,
              chromosome = sequencing.chromosome
      }
  }

  File filtered_gatk_vcf = select_first([ApplyRandomFilters.gatk_vcf_filt,  GatkGenotyping.vcf_norm])
  File filtered_gatk_vcf_bamcounts = select_first([ApplyRandomFilters.gatk_vcf_bam_counts_filt, GatkGenotyping.vcf_norm_bamcounts])
  File filtered_freebayes_vcf = select_first([ApplyRandomFilters.freebayes_vcf_filt, FreebayesGenotyping.vcf_norm])
  File filtered_freebayes_vcf_bamcounts = select_first([ApplyRandomFilters.freebayes_vcf_bam_counts_filt, FreebayesGenotyping.vcf_norm_bamcounts])
  
  call utils.GetMarkersPos {
    input:
      true_vcf = CreateAlignmentFromSimulation.true_vcf,
      filtered_gatk_vcf = filtered_gatk_vcf,
      filtered_gatk_vcf_bamcounts = filtered_gatk_vcf_bamcounts,
      filtered_freebayes_vcf = filtered_freebayes_vcf,
      filtered_freebayes_vcf_bamcounts = filtered_freebayes_vcf_bamcounts,
      depth = sequencing.depth,
      seed = family.seed
  }

  PopulationAnalysis gatk_processing = {"method": "gatk", "vcf": filtered_gatk_vcf, "bam": filtered_gatk_vcf_bamcounts}
  PopulationAnalysis freebayes_processing = {"method": "freebayes", "vcf": filtered_freebayes_vcf, "bam": filtered_freebayes_vcf_bamcounts}

  scatter (analysis in [gatk_processing, freebayes_processing]){

    Map[String, File] vcfs = {"vcf": analysis.vcf, "bam": analysis.bam}

    scatter (origin in ["vcf", "bam"]){

        call utils.SplitMarkers as splitgeno{
             input:
                vcf_file = vcfs[origin]
        }

        call genotyping.onemapMaps as updogMaps {
          input:
            simu_onemap_obj = truth_vcf.onemap_obj,
            vcf_file = splitgeno.biallelics,
            genotyping_program = "updog",
            ref_alt_alleles = CreateAlignmentFromSimulation.ref_alt_alleles,
            simulated_phases = CreateAlignmentFromSimulation.simulated_phases,
            SNPCall_program = analysis.method,
            CountsFrom = origin,
            cross = family.cross,
            max_cores = max_cores,
            simu_vcfR = truth_vcf.vcfR_obj,
            seed = family.seed,
            depth = sequencing.depth,
            multiallelics = sequencing.multiallelics,
            multiallelics_file = splitgeno.multiallelics
        }

        call genotyping.onemapMaps as supermassaMaps {
          input:
            simu_onemap_obj = truth_vcf.onemap_obj,
            vcf_file = splitgeno.biallelics,
            genotyping_program = "supermassa",
            ref_alt_alleles = CreateAlignmentFromSimulation.ref_alt_alleles,
            simulated_phases = CreateAlignmentFromSimulation.simulated_phases,
            SNPCall_program = analysis.method,
            CountsFrom = origin,
            cross = family.cross,
            max_cores = max_cores,
            simu_vcfR = truth_vcf.vcfR_obj,
            seed = family.seed,
            depth = sequencing.depth,
            multiallelics = sequencing.multiallelics,
            multiallelics_file = splitgeno.multiallelics
        }

        call genotyping.onemapMaps as polyradMaps {
          input:
            simu_onemap_obj = truth_vcf.onemap_obj,
            vcf_file = splitgeno.biallelics,
            genotyping_program = "polyrad",
            ref_alt_alleles = CreateAlignmentFromSimulation.ref_alt_alleles,
            simulated_phases = CreateAlignmentFromSimulation.simulated_phases,
            SNPCall_program = analysis.method,
            CountsFrom = origin,
            cross = family.cross,
            max_cores = max_cores,
            simu_vcfR = truth_vcf.vcfR_obj,
            seed = family.seed,
            depth = sequencing.depth,
            multiallelics = sequencing.multiallelics,
            multiallelics_file = splitgeno.multiallelics
        }
    }

    call utils.SplitMarkers as splitvcf{
         input:
           vcf_file = analysis.vcf
    }

    call utils.SplitMarkers as splitbam{
         input:
           vcf_file = analysis.bam
    }

    call gusmap.gusmapMaps {
      input:
        simu_onemap_obj = truth_vcf.onemap_obj,
        vcf_file = splitvcf.biallelics,
        new_vcf_file = splitbam.biallelics,
        SNPCall_program = analysis.method,
        GenotypeCall_program = "gusmap",
        ref_alt_alleles = CreateAlignmentFromSimulation.ref_alt_alleles,
        simulated_phases = CreateAlignmentFromSimulation.simulated_phases,
        seed = family.seed,
        depth = sequencing.depth,
        max_cores = max_cores
    }

    call snpcaller.SNPCallerMaps{
      input:
        simu_onemap_obj = truth_vcf.onemap_obj,
        vcf_file = splitvcf.biallelics,
        ref_alt_alleles = CreateAlignmentFromSimulation.ref_alt_alleles,
        simulated_phases = CreateAlignmentFromSimulation.simulated_phases,
        cross = family.cross,
        SNPCall_program = analysis.method,
        GenotypeCall_program = "SNPCaller",
        CountsFrom = "vcf",
        simu_vcfR = truth_vcf.vcfR_obj,
        seed = family.seed,
        depth = sequencing.depth,
        max_cores = max_cores,
        multiallelics = sequencing.multiallelics,
        multiallelics_file = splitvcf.multiallelics
    }
  }

  # Compress files
  call JointReports {
    input:
      SNPCaller                 = SNPCallerMaps.tar_gz_report,
      updog                     = flatten(updogMaps.tar_gz_report),
      polyrad                   = flatten(polyradMaps.tar_gz_report),
      supermassa                = flatten(supermassaMaps.tar_gz_report),
      gusmap                    = gusmapMaps.tar_gz_report,
      GATK_eval                 = GatkGenotyping.vcfEval,
      Freebayes_eval            = FreebayesGenotyping.vcfEval,
      multiallelics_file        = splitvcf.multiallelics,
      max_cores                 = max_cores,
      seed                      = family.seed,
      depth                     = sequencing.depth
  }

  output {
    File data1_depths_geno_prob   = JointReports.data1_depths_geno_prob
    File data2_maps               = JointReports.data2_maps
    File data3_filters            = JointReports.data3_filters
    File data4_times              = JointReports.data4_times
    File data5_SNPCall_efficiency = JointReports.data5_SNPCall_efficiency
    File data6_RDatas             = JointReports.data6_RDatas
    File data7_gusmap             = JointReports.data7_gusmap
    File data8_names              = JointReports.data8_names
    File data10_counts            = JointReports.data10_counts
    File simu_haplo               = CreateAlignmentFromSimulation.simu_haplo
    File Plots                    = GatkGenotyping.Plots
    File positions                = GetMarkersPos.positions
  }
}

task JointReports{
  input {
    Array[File] SNPCaller                
    Array[File] updog                     
    Array[File] polyrad                   
    Array[File] supermassa         
    Array[File] gusmap
    Array[File] multiallelics_file                    
    File Freebayes_eval
    File GATK_eval
    Int max_cores
    Int seed
    Int depth
  }

  command <<<
     R --vanilla --no-save <<RSCRIPT

      # packages
      library(tidyr)
      library(stringr)
      library(vroom)
      library(largeList)
      library(vcfR)

      SNPCaller  <- str_split("~{sep=";" SNPCaller}", ";", simplify = T)
      updog      <- str_split("~{sep=";" updog}", ";", simplify = T)
      polyrad    <- str_split("~{sep=";" polyrad}", ";", simplify = T)
      supermassa <- str_split("~{sep=";" supermassa}", ";", simplify = T)
      gusmap <- str_split("~{sep=";" gusmap}", ";", simplify = T)

      files <- list(SNPCaller, updog, polyrad, supermassa, gusmap)

      direc <- c("/maps/", "/filters/", "/errors/", "/times/", "/RDatas/")

      path_dir <- tempdir()
      system(paste0("mkdir ", paste0(path_dir, direc, collapse = " ")))
      for(i in 1:length(files)){
        for(j in 1:length(files[[i]])){
            untar(files[[i]][[j]], exdir = path_dir)
            list_files <- untar(files[[i]][[j]], exdir = path_dir, list = T)
            system(paste0("mv ",path_dir, "/",list_files[1], "*_map_report.tsv.gz ", path_dir, "/maps"))
            system(paste0("mv ",path_dir, "/",list_files[1], "*_times_report.tsv.gz ", path_dir, "/times"))
            system(paste0("mv ",path_dir, "/",list_files[1], "*.RData ", path_dir, "/RDatas"))
            if(!grepl("gusmap", list_files[1])){
              system(paste0("mv ",path_dir, "/",list_files[1], "*_filters_report.tsv.gz ", path_dir, "/filters"))
              system(paste0("mv ",path_dir, "/",list_files[1], "*_errors_report.tsv.gz ", path_dir, "/errors"))
            }
        }
      }

      direc_tsv <- direc[-5]
      tsvs <- list()
      for(i in 1:length(direc_tsv)){
        files <- system(paste0("ls ", path_dir, direc_tsv[i]), intern = T)
        files <- paste0(path_dir, direc_tsv[i], files)
        tsvs[[i]] <- vroom(files, num_threads = ~{max_cores})
      }

      # Add multiallelics tag
      vcf_multi <- str_split("~{sep = ";" multiallelics_file}", ";", simplify = T)

      vcfs <- lapply(as.list(vcf_multi), read.vcfR)
      multi_names_seed <- lapply(vcfs, function(x) paste0(x@fix[,1], "_", x@fix[,2]))

      snpcall_names <- rep(NA, length(vcfs))
      snpcall_names[which(sapply(vcfs, function(x) any(grep("gatk",x@meta))))] <- "gatk"
      snpcall_names[which(sapply(vcfs, function(x) any(grep("freebayes",x@meta))))] <- "freebayes"

      for(i in 1:length(snpcall_names)){
        tsvs[[1]][,"real.mks"][which(tsvs[[1]][,"SNPCall"] == snpcall_names[i] &  
                                      tsvs[[1]][,"mk.name"] %in% multi_names_seed[[i]])] <- "multiallelic"
      }

      library(gsalib)
      df <- gsa.read.gatkreport("~{Freebayes_eval}")
      eval1 <- cbind(SNPCall = "Freebayes", seed = 500, depth = 20, df[["CompOverlap"]])
      count1 <- cbind(SNPCall = "Freebayes", seed = 500, depth = 20, df[["CountVariants"]])

      df <- gsa.read.gatkreport("~{GATK_eval}")
      eval2 <- cbind(SNPCall = "GATK", seed = 500, depth = 20, df[["CompOverlap"]])
      count2 <- cbind(SNPCall = "GATK", seed = 500, depth = 20, df[["CountVariants"]])

      df <- rbind(eval1, eval2)
      vroom_write(df, "data5_SNPCall_efficiency.tsv.gz", num_threads = ~{max_cores})

      df <- rbind(count1, count2)
      vroom_write(df, "data10_CountVariants.tsv.gz", num_threads = ~{max_cores})

      rdatas_files <- paste0(path_dir, "/RDatas/",list.files(paste0(path_dir, "/RDatas/")))

      all_RDatas <- list()
      for(i in 1:length(rdatas_files)){
        map_temp <- load(rdatas_files[i])
        all_RDatas[[i]] <- get(map_temp)
      }
      all_RDatas <- unlist(all_RDatas, recursive = F)

      gusmap_RDatas <- all_RDatas[grep("gusmap", names(all_RDatas))]
      RDatas <- all_RDatas[-grep("gusmap", names(all_RDatas))]

      #   # Converting onemap sequencig objects to list. LargeList do not accept other class
      #   # Also because of this gusmap is separated, because the developers worked with enviroments, not classes
      
      for(i in 1:length(RDatas)){
        class(RDatas[[i]]) <- "list"
      }
      
      saveList(RDatas, file = "data6_RDatas.llo", append=FALSE, compress=TRUE)
      
      new_names <- names(all_RDatas)
      vroom_write(as.data.frame(new_names), "names.tsv.gz")
      save(gusmap_RDatas, file = "gusmap_RDatas.RData")

      # Outputs
      vroom_write(tsvs[[3]], "data1_depths_geno_prob.tsv.gz", num_threads = ~{max_cores})
      vroom_write(tsvs[[1]], "data2_maps.tsv.gz", num_threads = ~{max_cores})
      vroom_write(tsvs[[2]], "data3_filters.tsv.gz", num_threads = ~{max_cores})
      vroom_write(tsvs[[4]], "data4_times.tsv.gz", num_threads = ~{max_cores})

     RSCRIPT
  >>>

  runtime {
    docker:"cristaniguti/reads2map:0.0.1"
    # preemptible: 3
    # memory: "3 GB"
    # cpu: 4
    job_name:"JointReports"
    node:"--nodes=1"
    mem:"--mem=15GB"
    tasks:"--ntasks-per-node=11"
    time:"01:40:00"
    maxRetries: 5
  }

  output {
    File data1_depths_geno_prob = "data1_depths_geno_prob.tsv.gz"
    File data2_maps = "data2_maps.tsv.gz"
    File data3_filters = "data3_filters.tsv.gz"
    File data4_times   = "data4_times.tsv.gz"
    File data5_SNPCall_efficiency = "data5_SNPCall_efficiency.tsv.gz"
    File data6_RDatas  = "data6_RDatas.llo"
    File data7_gusmap  = "gusmap_RDatas.RData"
    File data8_names   = "names.tsv.gz"
    File data10_counts  = "data10_CountVariants.tsv.gz"
  }
}