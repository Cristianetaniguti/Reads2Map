version 1.0

import "structs/empiricalS.wdl"
import "tasks/create_alignment_from_families_files.wdl" as fam
import "tasks/gatk_genotyping.wdl" as gatk
import "tasks/freebayes_genotyping.wdl" as freebayes
import "tasks/utils.wdl" as utils
import "tasks/utilsR.wdl" as utilsR
import "tasks/default_maps_emp.wdl" as default
import "tasks/snpcaller_maps_emp.wdl" as snpcaller
import "tasks/updog_maps_emp.wdl" as updog
import "tasks/polyrad_maps_emp.wdl" as polyrad
import "tasks/supermassa_maps_emp.wdl" as supermassa
import "tasks/gusmap_maps_emp.wdl" as gusmap

workflow EmpiricalReads {

  input {
    Dataset dataset
    ReferenceFasta references
  }

  call fam.CreateAlignmentFromFamilies {
    input:
      families_info=dataset.samples_info,
      references=references
  }

  call gatk.GatkGenotyping {
    input:
      alignments=CreateAlignmentFromFamilies.alignments,
      references=references,
      program="gatk"
  }

  call freebayes.FreebayesGenotyping {
    input:
      alignments=CreateAlignmentFromFamilies.alignments,
      bam=CreateAlignmentFromFamilies.bam,
      bai=CreateAlignmentFromFamilies.bai,
      references=references,
      program="freebayes"
  }

  call utils.BamCounts4Onemap {
    input:
      sampleName=CreateAlignmentFromFamilies.names,
      freebayes_counts=FreebayesGenotyping.counts,
      gatk_counts=GatkGenotyping.counts
  }

  Array[String] methods = ["gatk", "freebayes"]
  Array[File] vcfs = [GatkGenotyping.vcf, FreebayesGenotyping.vcf]
  Array[Pair[String, File]] program_and_vcf = zip(methods, vcfs)
  
  scatter (vcf in program_and_vcf){
    call utilsR.vcf2onemap{
      input:
        vcf_file = vcf.right,
        cross = dataset.cross,
        SNPCall_program = vcf.left,
        parent1 = dataset.parent1,
        parent2 = dataset.parent2
    }
    
    call default.DefaultMaps{
      input:
        onemap_obj = vcf2onemap.onemap_obj,
        vcfR_obj = vcf2onemap.vcfR_obj,
        parent1 = dataset.parent1,
        parent2 = dataset.parent2,
        SNPCall_program = vcf.left,
        GenotypeCall_program = "default",
        CountsFrom = "vcf",
        cross = dataset.cross,
        chromosome = dataset.chromosome
    }
    
    call snpcaller.SNPCallerMaps{
      input:
        onemap_obj = vcf2onemap.onemap_obj,
        vcf_file = vcf.right,
        cross = dataset.cross,
        SNPCall_program = vcf.left,
        GenotypeCall_program = "SNPCaller",
        CountsFrom = "vcf",
        parent1 = dataset.parent1,
        parent2 = dataset.parent2,
        chromosome = dataset.chromosome
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
            onemap_obj = vcf2onemap.onemap_obj,
            vcf_file = vcf_counts.right, 
            SNPCall_program = vcf.left,
            GenotypeCall_program = "updog",
            CountsFrom = vcf_counts.left,
            cross = dataset.cross,
            parent1 = dataset.parent1,
            parent2 = dataset.parent2,
            chromosome = dataset.chromosome
        }
        
        call supermassa.SupermassaMaps{
          input:
            onemap_obj = vcf2onemap.onemap_obj,            
            vcf_file = vcf_counts.right, 
            SNPCall_program = vcf.left,
            GenotypeCall_program = "supermassa",
            CountsFrom = vcf_counts.left,
            cross = dataset.cross,
            parent1 = dataset.parent1,
            parent2 = dataset.parent2,
            chromosome = dataset.chromosome
        }
        
        call polyrad.PolyradMaps{
          input:
            onemap_obj = vcf2onemap.onemap_obj,            
            vcf_file = vcf_counts.right, 
            SNPCall_program = vcf.left,
            GenotypeCall_program = "polyrad",
            CountsFrom = vcf_counts.left,
            cross = dataset.cross,
            parent1 = dataset.parent1,
            parent2 = dataset.parent2,
            chromosome = dataset.chromosome
        }
      }
      call gusmap.GusmapMaps{
        input:
          vcf_file = vcf.right,
          new_vcf_file = BamDepths2Vcf.bam_vcf,
          SNPCall_program = vcf.left,
          GenotypeCall_program = "gusmap",
          parent1 = dataset.parent1,
          parent2 = dataset.parent2,
          chromosome = dataset.chromosome
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
      
     library(tidyr)
     library(largeList)
      
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

      gusmap_RDatas <- all_RDatas[grep("gusmap", names(all_RDatas))]
      RDatas <- all_RDatas[-grep("gusmap", names(all_RDatas))]
      
      # Converting OneMap sequencig objects to list. LargeList do not accept other class
      # Also because of this gusmap is separated, because the developers worked with enviroments, not classes
      
      for(i in 1:length(RDatas)){
        class(RDatas[[i]]) <- "list"
      }
      
      saveList(RDatas, file = "onemap_RDatas.llo", append=FALSE, compress=TRUE)
      
      new_names <- names(all_RDatas)
      saveRDS(new_names, file = "names.rds")
      save(gusmap_RDatas, file = "gusmap_RDatas.RData")
          
     RSCRIPT
  >>>
  runtime{
    docker:"taniguti/onemap"
    time:"05:00:00"
    mem:"--nodes=1"
    cpu:1
  }
  
  output{
    File onemap_RDatas = "onemap_RDatas.RData"
    File gusmap_RDatas = "gusmap_RDatas.RData"
    File all_maps = "all_maps.txt"
    File all_filters = "all_filters.txt"
    File all_errors = "all_errors.txt"
    File all_times = "all_times.txt"
  }
}

