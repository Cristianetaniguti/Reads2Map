version 1.0

import "../structs/reads_simuS.wdl"
import "./create_alignment_from_read_simulations.wdl" as simulation
import "./gatk_genotyping.wdl" as gatk
import "./freebayes_genotyping.wdl" as freebayes
import "./utils.wdl" as utils
import "./utilsR.wdl" as utilsR
import "./simulated_map.wdl" as simulated_map
import "./default_maps.wdl" as default
import "./snpcaller_maps.wdl" as snpcaller
import "./updog_maps.wdl" as updog
import "./polyrad_maps.wdl" as polyrad
import "./supermassa_maps.wdl" as supermassa
import "./gusmap_maps.wdl" as gusmap


workflow reads_simu {

  input {
    ReferenceFasta references
    Family family
    Profiles profiles
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
      program="gatk"
  }

  call freebayes.FreebayesGenotyping {
    input:
      alignments=CreateAlignmentFromSimulation.alignments,
      bam=CreateAlignmentFromSimulation.bam,
      bai=CreateAlignmentFromSimulation.bai,
      references=references,
      program="freebayes"
  }

  call utils.CalculateVcfMetrics {
    input:
      freebayesVCF  = FreebayesGenotyping.vcf,
      gatkVCF       = GatkGenotyping.vcf,
      tot_mks       = CreateAlignmentFromSimulation.total_markers,
      maternal_trim = CreateAlignmentFromSimulation.maternal_trim,
      seed          = family.seed,
      depth         = family.depth
  }

  call utils.BamCounts4Onemap {
    input:
      sampleName       = CreateAlignmentFromSimulation.names,
      freebayes_counts = FreebayesGenotyping.counts,
      gatk_counts      = GatkGenotyping.counts
  }
  
  call simulated_map.SimulatedMap{
    input:
      vcf_simu = CreateAlignmentFromSimulation.true_vcf,
      cross = family.cross
  }
  
  Array[String] methods                     = ["gatk", "freebayes"]
  Array[File] vcfs                          = [GatkGenotyping.vcf, FreebayesGenotyping.vcf]
  Array[Pair[String, File]] program_and_vcf = zip(methods, vcfs)

  scatter (vcf in program_and_vcf){
  
    call utilsR.vcf2onemap{
      input:
        vcf_file = vcf.right,
        cross = family.cross,
        SNPCall_program = vcf.left
    }
  
    call default.DefaultMaps{
      input:
        simu_onemap_obj = SimulatedMap.simu_onemap_obj,
        onemap_obj = vcf2onemap.onemap_obj,
        tot_mks = CreateAlignmentFromSimulation.total_markers,
        real_phases = CreateAlignmentFromSimulation.real_phases,
        SNPCall_program = vcf.left,
        GenotypeCall_program = "default",
        CountsFrom = "vcf",
        cMbyMb = family.cmBymb,
        cross = family.cross
    }
    
    call snpcaller.SNPCallerMaps{
      input:
        simu_onemap_obj = SimulatedMap.simu_onemap_obj,
        onemap_obj = vcf2onemap.onemap_obj,
        vcf_file = vcf.right,
        tot_mks = CreateAlignmentFromSimulation.total_markers,
        real_phases = CreateAlignmentFromSimulation.real_phases,
        cross = family.cross,
        SNPCall_program = vcf.left,
        GenotypeCall_program = "SNPCaller",
        CountsFrom = "vcf",
        cMbyMb = family.cmBymb
    }
    
    call utilsR.BamDepths2Vcf{
      input:
        SNPCall_program = vcf.left,
        vcf_file = vcf.right,
        freebayes_ref_bam = BamCounts4Onemap.freebayes_ref_bam,
        freebayes_alt_bam = BamCounts4Onemap.freebayes_alt_bam,
        gatk_ref_bam = BamCounts4Onemap.gatk_ref_bam,
        gatk_alt_bam = BamCounts4Onemap.gatk_alt_bam,
        gatk_example_alleles      = BamCounts4Onemap.gatk_example_alleles,
        freebayes_example_alleles = BamCounts4Onemap.freebayes_example_alleles
    }
    
    Array[String] counts     = ["vcf", "bam"]
    Array[File] vcfs_counts  = [vcf.right, BamDepths2Vcf.bam_vcf]
    Array[Pair[String, File]] counts_and_vcf = zip(counts, vcfs_counts)
    
    scatter(vcf_counts in counts_and_vcf){
        call updog.UpdogMaps{
          input:
            simu_onemap_obj = SimulatedMap.simu_onemap_obj,
            onemap_obj = vcf2onemap.onemap_obj,
            vcf_file = vcf_counts.right, 
            tot_mks = CreateAlignmentFromSimulation.total_markers,
            real_phases = CreateAlignmentFromSimulation.real_phases,
            SNPCall_program = vcf.left,
            GenotypeCall_program = "updog",
            CountsFrom = vcf_counts.left,
            cMbyMb = family.cmBymb,
            cross = family.cross
        }
        
        call supermassa.SupermassaMaps{
          input:
            simu_onemap_obj = SimulatedMap.simu_onemap_obj,
            onemap_obj = vcf2onemap.onemap_obj,            
            vcf_file = vcf_counts.right, 
            tot_mks = CreateAlignmentFromSimulation.total_markers,
            real_phases = CreateAlignmentFromSimulation.real_phases,
            SNPCall_program = vcf.left,
            GenotypeCall_program = "supermassa",
            CountsFrom = vcf_counts.left,
            cMbyMb = family.cmBymb,
            cross = family.cross
        }
        
        call polyrad.PolyradMaps{
          input:
            simu_onemap_obj = SimulatedMap.simu_onemap_obj,
            onemap_obj = vcf2onemap.onemap_obj,            
            vcf_file = vcf_counts.right, 
            tot_mks = CreateAlignmentFromSimulation.total_markers,
            real_phases = CreateAlignmentFromSimulation.real_phases,
            SNPCall_program = vcf.left,
            GenotypeCall_program = "polyrad",
            CountsFrom = vcf_counts.left,
            cMbyMb = family.cmBymb,
            cross = family.cross
        }
      }
      call gusmap.GusmapMaps{
        input:
          simu_onemap_obj = SimulatedMap.simu_onemap_obj,
          vcf_file = vcf.right,
          new_vcf_file = BamDepths2Vcf.bam_vcf,
          SNPCall_program = vcf.left,
          GenotypeCall_program = "gusmap",
          tot_mks = CreateAlignmentFromSimulation.total_markers,
          real_phases = CreateAlignmentFromSimulation.real_phases,
          cMbyMb = family.cmBymb
      }
  }
  call JointReports{
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
  
  call CreateTables{
  input:
    depth                     = family.depth,
    seed                      = family.seed,
    gatk_ref_depth            = CalculateVcfMetrics.gatk_ref_depth,
    gatk_ref_depth_bam        = BamCounts4Onemap.gatk_ref_bam,
    gatk_alt_depth            = CalculateVcfMetrics.gatk_alt_depth,
    gatk_alt_depth_bam        = BamCounts4Onemap.gatk_alt_bam,
    freebayes_ref_depth_bam   = BamCounts4Onemap.freebayes_ref_bam,
    freebayes_alt_depth_bam   = BamCounts4Onemap.freebayes_alt_bam,
    freebayes_ref_depth       = CalculateVcfMetrics.freebayes_ref_depth,
    freebayes_alt_depth       = CalculateVcfMetrics.freebayes_alt_depth,
    all_maps                  = JointReports.all_maps,
    all_errors                = JointReports.all_errors,
    all_filters               = JointReports.all_filters,
    all_times                 = JointReports.all_times,
    all_RDatas                = JointReports.all_RDatas
  }
  
  output{
    File data1_depths_geno_prob   = CreateTables.data1_depths_geno_prob
    File data2_maps               = CreateTables.data2_maps
    File data3_filters            = CreateTables.data3_filters
    File data5_SNPcall_efficiency = CalculateVcfMetrics.data5_SNPcall_efficiency
    File data4_times              = CreateTables.data4_times
    File data6_RDatas             = CreateTables.data6_RDatas
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
      system("cat ~{sep= ' ' default_maps_report} ~{sep= ' ' SNPCaller_maps_report} ~{sep= ' ' Updog_maps_report} ~{sep= ' ' Polyrad_maps_report} ~{sep= ' ' Supermassa_maps_report}  ~{sep= ' ' Gusmap_maps_report} > all_maps.txt")
     
      system("cat ~{sep= ' ' default_filters_report} ~{sep= ' ' SNPCaller_filters_report} ~{sep= ' ' Updog_filters_report} ~{sep= ' ' Polyrad_filters_report} ~{sep= ' ' Supermassa_filters_report}  > all_filters.txt")
     
      system("cat ~{sep= ' '  default_errors_report} ~{sep= ' ' SNPCaller_errors_report} ~{sep= ' ' Updog_errors_report} ~{sep= ' ' Polyrad_errors_report} ~{sep= ' ' Supermassa_errors_report} > all_errors.txt")
      
      system("cat ~{sep= ' '  default_times} ~{sep= ' ' SNPCaller_times} ~{sep= ' ' Updog_times} ~{sep= ' ' Polyrad_times} ~{sep= ' ' Supermassa_times}  ~{sep= ' ' Gusmap_times} > all_times.txt")
      
      system("cp ~{sep= ' ' default_RDatas} ~{sep= ' ' SNPCaller_RDatas}  ~{sep= ' ' Updog_RDatas}  ~{sep= ' ' Polyrad_RDatas}  ~{sep= ' ' Supermassa_RDatas} ~{sep= ' ' Gusmap_RDatas} .")
      
     Genocall <- c("default", "SNPCaller", "updog", "supermassa", "polyrad", "gusmap",
                   "default0.05", "updog0.05", "supermassa0.05", "polyrad0.05")
     fake <- c(TRUE, FALSE)
     CountsFrom <- c("vcf", "bam")
     SNPCall <- c("gatk", "freebayes")
     
      all_RDatas <- list()
      z <- 1
      names_RDatas <- vector()
      for(i in 1:length(Genocall)){
        for(j in 1:length(CountsFrom)){
          for(w in 1:length(fake)){
            for(y in 1:length(SNPCall)){
              if(CountsFrom[j] == "bam" & (Genocall[i] == "default" | Genocall[i] == "SNPCaller" |  Genocall[i] == "default0.05")){
            } else {
              map_temp <- load(paste0("map_",SNPCall[y],"_",CountsFrom[j], "_", Genocall[i],".RData"))
              all_RDatas[[z]] <- get(map_temp)
              names_RDatas <- c(names_RDatas, paste0("map_",SNPCall[y],"_",CountsFrom[j], "_", Genocall[i]))
              z <- z+1
              }
            }
          }
        }
      }
      save(all_RDatas, file= "all_RDatas.RData")
          
     RSCRIPT
  >>>
  runtime{
    docker:"taniguti/onemap"
    time:"-5:00:00"
    mem:"--nodes=1"
    cpu:1
  }
  
  output{
    File all_RDatas = "all_RDatas.RData"
    File all_maps = "all_maps.txt"
    File all_filters = "all_filters.txt"
    File all_errors = "all_errors.txt"
    File all_times = "all_times.txt"
  }
}


task CreateTables{
  input{
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
    File all_maps
    File all_errors
    File all_filters
    File all_times
    File all_RDatas
  }
  
  command <<<
    
    R --vanilla --no-save <<RSCRIPT
  
      library(reshape2)
      library(vcfR)
      
      depth <- ~{depth}
      seed <- ~{seed}
      
      all_errors <- read.table("~{all_errors}")
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
      
      all_errors <- cbind(seed, depth, data1)
      
      ########################################################
      # Table2: seed; CountsFrom; ErrorProb; SNPcall; MK; rf; phases; real_phases
      ########################################################
      # Add seed and mean depth information
      all_maps <- read.table("~{all_maps}")
      all_maps <- cbind(seed=seed, depth=depth, all_maps)
    
      ##########################################################################
      # Table3: CountsFrom; seed; SNPcall; GenoCall; n_mks; distorted; redundant
      ##########################################################################
      all_filters <- read.table("~{all_filters}")
      all_filters <- cbind(seed=seed, depth=depth, all_filters)
    
      ###########################################################################
      # Table4: CountsFrom; seed; SNPcall; GenoCall
      ###########################################################################
      all_times <- read.table("~{all_times}", stringsAsFactors = F)
      
      temp <- strsplit(all_times[,2], "_")
      temp <- do.call(rbind, temp)
      temp <- temp[,-1]
      colnames(temp) <- c("SNPCall", "CountsFrom", "GenoCall", "fake")
      all_times <- cbind(seed=seed, depth=depth, temp, time= all_times[,3])
    
      ###########################################################################
      # Table6: list of RDatas with name CountsFrom; seed; SNPcall; GenoCall
      ###########################################################################
      load("~{all_RDatas}")
      RDatas <- unlist(all_RDatas, recursive = F)
      all_names <- unlist(lapply(all_RDatas, names))
      new_names <- paste0(seed, "_", depth, "_", all_names)
      names(RDatas) <- new_names
    
      saveRDS(all_errors, file = "data1_depths_geno_prob.rds")
      saveRDS(all_maps, file = "data2_maps.rds")
      saveRDS(all_filters, file = "data3_filters.rds")
      saveRDS(all_times, file= "data4_times.rds")
      save(RDatas, file = "data6_RDatas.RData")
  
  RSCRIPT
  >>>
    
    runtime{
      docker:"taniguti/onemap"
      time:"05:00:00"
      mem:"--nodes=1"
      cpu:1
    }
  
  output{
    File data1_depths_geno_prob = "data1_depths_geno_prob.rds"
    File data2_maps = "data2_maps.rds"
    File data3_filters = "data3_filters.rds"
    File data4_times   = "data4_times.rds"
    File data6_RDatas  = "data6_RDatas.RData"
  }
}


