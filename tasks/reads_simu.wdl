version 1.0


import "./create_alignment_from_read_simulations.wdl" as simulation
import "./gatk_genotyping.wdl" as gatk
import "./freebayes_genotyping.wdl" as freebayes
import "./utils.wdl" as utils
import "./utilsR.wdl" as utilsR
import "./simulated_map.wdl" as simulated_map
import "./default_maps.wdl" as default
import "./snpcaller_maps.wdl" as snpcaller
import "./gusmap_maps.wdl" as gusmap
import "./genotyping-simulated.wdl" as genotyping


struct PopulationAnalysis {
    String method
    File vcf
    File bam
    File? multi
}

workflow reads_simu {

  input {
    Reference references
    Family family
    Profiles profiles
    SplitVCF splitvcf
    String? filters
  }

  call simulation.CreateAlignmentFromSimulation {
    input:
      references=references,
      family=family,
      profiles=profiles
  }

  call gatk.GatkGenotyping {
    input:
      alignments=CreateAlignmentFromSimulation.alignments,
      references=references,
      program="gatk",
      splitvcf = splitvcf,
      sampleNames = CreateAlignmentFromSimulation.names
  }

  call freebayes.FreebayesGenotyping {
    input:
      alignments=CreateAlignmentFromSimulation.alignments,
      bam=CreateAlignmentFromSimulation.bam,
      bai=CreateAlignmentFromSimulation.bai,
      references=references,
      program="freebayes",
      splitvcf = splitvcf,
      sampleNames = CreateAlignmentFromSimulation.names
  }

  call utils.CalculateVcfMetrics {
    input:
      freebayesVCF  = FreebayesGenotyping.vcf_bi,
      gatkVCF       = GatkGenotyping.vcf_bi,
      tot_mks       = CreateAlignmentFromSimulation.total_markers,
      maternal_trim = CreateAlignmentFromSimulation.maternal_trim,
      seed          = family.seed,
      depth         = family.depth
  }

  call simulated_map.SimulatedMap{
    input:
      vcf_simu = CreateAlignmentFromSimulation.true_vcf,
      cross = family.cross
  }

    if (defined(filters)) {
        call utils.ApplyRandomFilters {
            input:
                gatk_vcf = GatkGenotyping.vcf_bi,
                freebayes_vcf = FreebayesGenotyping.vcf_bi,
                gatk_vcf_bam_counts = GatkGenotyping.vcf_bi_bam_counts,
                freebayes_vcf_bam_counts = FreebayesGenotyping.vcf_bi_bam_counts,
                filters = filters,
                chromosome = splitvcf.chromosome
        }
    }

    File filtered_gatk_vcf = select_first([ApplyRandomFilters.gatk_vcf_filt,  GatkGenotyping.vcf_bi])
    File filtered_gatk_vcf_bamcounts = select_first([ApplyRandomFilters.gatk_vcf_bam_counts_filt, GatkGenotyping.vcf_bi_bam_counts])
    File filtered_freebayes_vcf = select_first([ApplyRandomFilters.freebayes_vcf_filt, FreebayesGenotyping.vcf_bi])
    File filtered_freebayes_vcf_bamcounts = select_first([ApplyRandomFilters.freebayes_vcf_bam_counts_filt, FreebayesGenotyping.vcf_bi_bam_counts])


    PopulationAnalysis gatk_processing = {"multi":GatkGenotyping.vcf_multi, "method": "gatk", "vcf": filtered_gatk_vcf, "bam": filtered_gatk_vcf_bamcounts}
    PopulationAnalysis freebayes_processing = {"multi":FreebayesGenotyping.vcf_multi, "method": "freebayes", "vcf": filtered_freebayes_vcf, "bam": filtered_freebayes_vcf_bamcounts}


  scatter (analysis in [gatk_processing, freebayes_processing]){

    call utilsR.vcf2onemap{
      input:
        vcf_file = analysis.vcf,
        cross = family.cross,
        SNPCall_program = analysis.method,
        parent1 = "P1",
        parent2 = "P2"
    }

    if (family.multiallelics == "true"){
      call utilsR.MultiVcf2onemap{
         input:
            multi = analysis.multi,
            cross = family.cross,
            SNPCall_program = analysis.method,
            parent1 = "P1",
            parent2 = "P2",
      }
    }

    call default.DefaultMaps {
      input:
        onemap_obj = vcf2onemap.onemap_obj,
        simu_onemap_obj = SimulatedMap.simu_onemap_obj,
        tot_mks = CreateAlignmentFromSimulation.total_markers,
        real_phases = CreateAlignmentFromSimulation.real_phases,
        SNPCall_program = analysis.method,
        CountsFrom = "vcf",
        cMbyMb = family.cmBymb,
        multi_obj = MultiVcf2onemap.onemap_obj
    }

    call snpcaller.SNPCallerMaps{
      input:
        simu_onemap_obj = SimulatedMap.simu_onemap_obj,
        onemap_obj = vcf2onemap.onemap_obj,
        vcf_file = analysis.vcf,
        tot_mks = CreateAlignmentFromSimulation.total_markers,
        real_phases = CreateAlignmentFromSimulation.real_phases,
        cross = family.cross,
        SNPCall_program = analysis.method,
        GenotypeCall_program = "SNPCaller",
        CountsFrom = "vcf",
        cMbyMb = family.cmBymb,
        multi_obj = MultiVcf2onemap.onemap_obj
    }

    Map[String, File] vcfs = {"vcf": analysis.vcf, "bam": analysis.bam}

    scatter (origin in ["vcf", "bam"]){
        call genotyping.SnpBasedGenotypingSimulatedMaps as UpdogMaps {
          input:
            simu_onemap_obj = SimulatedMap.simu_onemap_obj,
            onemap_obj = vcf2onemap.onemap_obj,
            vcf_file = vcfs[origin],
            genotyping_program = "updog",
            tot_mks = CreateAlignmentFromSimulation.total_markers,
            real_phases = CreateAlignmentFromSimulation.real_phases,
            SNPCall_program = analysis.method,
            CountsFrom = origin,
            cMbyMb = family.cmBymb,
            cross = family.cross,
            multi_obj = MultiVcf2onemap.onemap_obj
        }

        call genotyping.SnpBasedGenotypingSimulatedMaps as SupermassaMaps {
          input:
            simu_onemap_obj = SimulatedMap.simu_onemap_obj,
            onemap_obj = vcf2onemap.onemap_obj,
            vcf_file = vcfs[origin],
            genotyping_program = "supermassa",
            tot_mks = CreateAlignmentFromSimulation.total_markers,
            real_phases = CreateAlignmentFromSimulation.real_phases,
            SNPCall_program = analysis.method,
            CountsFrom = origin,
            cMbyMb = family.cmBymb,
            cross = family.cross,
            multi_obj = MultiVcf2onemap.onemap_obj
        }

        call genotyping.SnpBasedGenotypingSimulatedMaps as PolyradMaps {
          input:
            simu_onemap_obj = SimulatedMap.simu_onemap_obj,
            onemap_obj = vcf2onemap.onemap_obj,
            vcf_file = vcfs[origin],
            genotyping_program = "polyrad",
            tot_mks = CreateAlignmentFromSimulation.total_markers,
            real_phases = CreateAlignmentFromSimulation.real_phases,
            SNPCall_program = analysis.method,
            CountsFrom = origin,
            cMbyMb = family.cmBymb,
            cross = family.cross,
            multi_obj = MultiVcf2onemap.onemap_obj
        }
      }

      call gusmap.GusmapMaps {
        input:
          simu_onemap_obj = SimulatedMap.simu_onemap_obj,
          vcf_file = analysis.vcf,
          new_vcf_file = analysis.bam,
          SNPCall_program = analysis.method,
          GenotypeCall_program = "gusmap",
          tot_mks = CreateAlignmentFromSimulation.total_markers,
          real_phases = CreateAlignmentFromSimulation.real_phases,
          cMbyMb = family.cmBymb
      }
  }

  call JointReports {
    input:
    default_RDatas            = flatten(DefaultMaps.RDatas),
    default_maps_report       = flatten(DefaultMaps.maps_report),
    default_filters_report    = flatten(DefaultMaps.filters_report),
    default_errors_report     = flatten(DefaultMaps.errors_report),
    default_times             = flatten(DefaultMaps.times),
    SNPCaller_RDatas          = SNPCallerMaps.RDatas,
    SNPCaller_maps_report     = SNPCallerMaps.maps_report,
    SNPCaller_filters_report  = SNPCallerMaps.filters_report,
    SNPCaller_errors_report   = SNPCallerMaps.errors_report,
    SNPCaller_times           = SNPCallerMaps.times,
    Updog_RDatas              = flatten(flatten(UpdogMaps.RDatas)),
    Updog_maps_report         = flatten(flatten(UpdogMaps.maps_report)),
    Updog_filters_report      = flatten(flatten(UpdogMaps.filters_report)),
    Updog_errors_report       = flatten(flatten(UpdogMaps.errors_report)),
    Updog_times               = flatten(flatten(UpdogMaps.times)),
    Polyrad_RDatas            = flatten(flatten(PolyradMaps.RDatas)),
    Polyrad_maps_report       = flatten(flatten(PolyradMaps.maps_report)),
    Polyrad_filters_report    = flatten(flatten(PolyradMaps.filters_report)),
    Polyrad_errors_report     = flatten(flatten(PolyradMaps.errors_report)),
    Polyrad_times             = flatten(flatten(PolyradMaps.times)),
    Supermassa_RDatas         = flatten(flatten(SupermassaMaps.RDatas)),
    Supermassa_maps_report    = flatten(flatten(SupermassaMaps.maps_report)),
    Supermassa_filters_report = flatten(flatten(SupermassaMaps.filters_report)),
    Supermassa_errors_report  = flatten(flatten(SupermassaMaps.errors_report)),
    Supermassa_times          = flatten(flatten(SupermassaMaps.times)),
    Gusmap_RDatas             = flatten(GusmapMaps.RDatas),
    Gusmap_maps_report        = flatten(GusmapMaps.maps_report),
    Gusmap_times              = flatten(GusmapMaps.times),
    depth                     = family.depth,
    seed                      = family.seed,
    gatk_ref_depth            = CalculateVcfMetrics.gatk_ref_depth,
    gatk_ref_depth_bam        = GatkGenotyping.ref_bam,
    gatk_alt_depth            = CalculateVcfMetrics.gatk_alt_depth,
    gatk_alt_depth_bam        = GatkGenotyping.alt_bam,
    freebayes_ref_depth_bam   = FreebayesGenotyping.ref_bam,
    freebayes_alt_depth_bam   = FreebayesGenotyping.alt_bam,
    freebayes_ref_depth       = CalculateVcfMetrics.freebayes_ref_depth,
    freebayes_alt_depth       = CalculateVcfMetrics.freebayes_alt_depth
  }

  output {
    File data1_depths_geno_prob   = JointReports.data1_depths_geno_prob
    File data2_maps               = JointReports.data2_maps
    File data3_filters            = JointReports.data3_filters
    File data5_SNPcall_efficiency = CalculateVcfMetrics.data5_SNPcall_efficiency
    File data4_times              = JointReports.data4_times
    File data6_RDatas             = JointReports.data6_RDatas
    File data7_gusmap             = JointReports.data7_gusmap
    File data8_names              = JointReports.data8_names
    File simu_haplo               = CreateAlignmentFromSimulation.simu_haplo
  }
}

task JointReports{
  input {
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
    Int depth
    Int seed
    File gatk_ref_depth
    File gatk_ref_depth_bam
    File gatk_alt_depth
    File gatk_alt_depth_bam
    File freebayes_ref_depth_bam
    File freebayes_alt_depth_bam
    File freebayes_ref_depth
    File freebayes_alt_depth
  }

  command <<<
     R --vanilla --no-save <<RSCRIPT

      # Joint reports
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
      library(reshape2)
      library(largeList)

      depth <- ~{depth}
      seed <- ~{seed}

      all_errors <- read.table("all_errors.txt")

      # Adding allele counts to the report
      gatk.ref.depth <- read.table("~{gatk_ref_depth}")

      chr <- unique(sapply(strsplit(rownames(gatk.ref.depth), "_"), "[",1))

      all_errors[,5] <- paste0(chr, "_",all_errors[,5])
      colnames(all_errors) <- c("SNPCall", "GenoCall", "CountsFrom", "ind", "mks", "gabGT", "methGT", "A", "AB", "BA", "B")

      files <- list(c("~{gatk_ref_depth}", "~{gatk_alt_depth}"), c("~{gatk_ref_depth_bam}", "~{gatk_alt_depth_bam}"),
                    c("~{freebayes_ref_depth}", "~{freebayes_alt_depth}"), c("~{freebayes_ref_depth_bam}",
                                                                               "~{freebayes_alt_depth_bam}"))
      ref_alt <- c("ref", "alt")
      SNPCall <- rep(c("gatk", "freebayes"), each=4)
      CountsFrom <- rep(c("vcf", "bam"),2, each=2)

      ########################################################################################
      # Table1: GenoCall; mks; ind; SNPcall; CountsFrom; alt; ref; gabGT; methGT; A; AB; BA; B
      ########################################################################################

      # Binding depth information for each genotype
      z <- 1
      data1 <- data.frame()
      for(j in 1:length(files)){
        for(i in 1:2){
          if(i == 1) df_meth <- all_errors
          alleles <- read.table(files[[j]][i])
          alleles <- cbind(mks=rownames(alleles), alleles)
          alleles <- melt(alleles, id.vars = c("mks"))
          colnames(alleles) <- c("mks", "ind", ref_alt[i])
          alleles <- cbind(SNPCall=SNPCall[z], CountsFrom = CountsFrom[z], alleles)
          z <- z+ 1
          df_meth <- merge(df_meth, alleles, by = c("SNPCall", "CountsFrom", "ind", "mks"))

        }
        data1 <- rbind(data1, df_meth)
      }

      # Adding seed and depth info
      all_errors <- cbind(seed, depth, data1)

      ########################################################
      # Table2: seed; CountsFrom; ErrorProb; SNPcall; MK; rf; phases; real_phases
      ########################################################
      # Add seed and mean depth information
      all_maps <- read.table("all_maps.txt")
      all_maps <- cbind(seed=seed, depth=depth, all_maps)

      ##########################################################################
      # Table3: CountsFrom; seed; SNPcall; GenoCall; n_mks; distorted; redundant
      ##########################################################################
      all_filters <- read.table("all_filters.txt")
      all_filters <- cbind(seed=seed, depth=depth, all_filters)

      ###########################################################################
      # Table4: CountsFrom; seed; SNPcall; GenoCall
      ###########################################################################
      all_times <- read.table("all_times.txt", stringsAsFactors = F)

      temp <- strsplit(all_times[,2], "_")
      temp <- do.call(rbind, temp)
      temp <- temp[,-1]
      colnames(temp) <- c("SNPCall", "CountsFrom", "GenoCall", "fake")
      all_times <- cbind(seed=seed, depth=depth, temp, time= all_times[,3])

      ###########################################################################
      # Table6: list of RDatas with name CountsFrom; seed; depth; SNPcall; GenoCall
      ###########################################################################
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
      all_RDatas <- unlist(all_RDatas, recursive = F)

      all_names <- names(all_RDatas)
      new_names <- paste0(seed, "_", depth, "_", all_names)
      names(all_RDatas) <- new_names

      # Outputs
      saveRDS(all_errors, file = "data1_depths_geno_prob.rds")
      saveRDS(all_maps, file = "data2_maps.rds")
      saveRDS(all_filters, file = "data3_filters.rds")
      saveRDS(all_times, file= "data4_times.rds")

      gusmap_RDatas <- all_RDatas[grep("gusmap", names(all_RDatas))]
      RDatas <- all_RDatas[-grep("gusmap", names(all_RDatas))]

      # Converting OneMap sequencig objects to list. LargeList do not accept other class
      # Also because of this gusmap is separated, because the developers worked with enviroments, not classes

      for(i in 1:length(RDatas)){
        class(RDatas[[i]]) <- "list"
      }

      saveList(RDatas, file = "data6_RDatas.llo", append=FALSE, compress=TRUE)

      # LargeList package limits to 16 letters each character, therefore, the entire names are stored in a separated vector
      saveRDS(new_names, file = "names.rds")
      save(gusmap_RDatas, file = "gusmap_RDatas.RData")

     RSCRIPT

  >>>

  runtime {
    docker:"cristaniguti/onemap_workflows"
    time:"48:00:00"
    mem:"--nodes=1"
    cpu:1
  }

  output {
    File data1_depths_geno_prob = "data1_depths_geno_prob.rds"
    File data2_maps = "data2_maps.rds"
    File data3_filters = "data3_filters.rds"
    File data4_times   = "data4_times.rds"
    File data6_RDatas  = "data6_RDatas.llo"
    File data7_gusmap  = "gusmap_RDatas.RData"
    File data8_names   = "names.rds"
  }
}
