version 1.0

import "../structs/reads_simuS.wdl"
import "./create_alignment_from_read_simulations.wdl" as simulation
import "./gatk_genotyping.wdl" as gatk
import "./freebayes_genotyping.wdl" as freebayes
import "./utils.wdl" as utils


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

  Array[String] methods                     = ["gatk", "freebayes"]
  Array[File] vcfs                          = [GatkGenotyping.vcf, FreebayesGenotyping.vcf]
  Array[Pair[String, File]] program_and_vcf = zip(methods, vcfs)

  scatter (vcf in program_and_vcf){
    call CreateMaps{
      input:
        tot_mks                   = CreateAlignmentFromSimulation.total_markers,
        simu_vcf                  = CreateAlignmentFromSimulation.true_vcf,
        methodName                = vcf.left,
        vcf_file                  = vcf.right,
        freebayes_ref_depth       = BamCounts4Onemap.freebayes_ref_bam,
        freebayes_alt_depth       = BamCounts4Onemap.freebayes_alt_bam,
        gatk_ref_depth            = BamCounts4Onemap.gatk_ref_bam,
        gatk_alt_depth            = BamCounts4Onemap.gatk_alt_bam,
        gatk_example_alleles      = BamCounts4Onemap.gatk_example_alleles,
        freebayes_example_alleles = BamCounts4Onemap.freebayes_example_alleles,
        cross                     = family.cross,
        real_phases               = CreateAlignmentFromSimulation.real_phases,
        cmBymb                    = family.cmBymb
    }
  }

  call CreateTables{
    input:
        depth                     = family.depth,
        seed                      = family.seed,
        tot_mks                   = CreateAlignmentFromSimulation.total_markers,
        gatk_ref_depth            = CalculateVcfMetrics.gatk_ref_depth,
        gatk_ref_depth_bam        = BamCounts4Onemap.gatk_ref_bam,
        gatk_alt_depth            = CalculateVcfMetrics.gatk_alt_depth,
        gatk_alt_depth_bam        = BamCounts4Onemap.gatk_alt_bam,
        freebayes_ref_depth_bam   = BamCounts4Onemap.freebayes_ref_bam,
        freebayes_alt_depth_bam   = BamCounts4Onemap.freebayes_alt_bam,
        freebayes_ref_depth       = CalculateVcfMetrics.freebayes_ref_depth,
        freebayes_alt_depth       = CalculateVcfMetrics.freebayes_alt_depth,
        all_maps                  = CreateMaps.all_maps,
        all_errors                = CreateMaps.all_errors,
        all_filters               = CreateMaps.all_filters,
        times                     = CreateMaps.times,
        all_RDatas                = CreateMaps.all_RDatas
    }

  output {
    File data1_depths_geno_prob   = CreateTables.data1_depths_geno_prob
    File data2_maps               = CreateTables.data2_maps
    File data3_filters            = CreateTables.data3_filters
    File data5_SNPcall_efficiency = CalculateVcfMetrics.data5_SNPcall_efficiency
    File data4_times              = CreateTables.data4_times
    File data6_RDatas             = CreateTables.data6_RDatas
  }
}

task CreateMaps {

 input {
    File tot_mks
    String methodName
    File simu_vcf
    File vcf_file
    File freebayes_ref_depth
    File freebayes_alt_depth
    File gatk_ref_depth
    File gatk_alt_depth
    File gatk_example_alleles
    File freebayes_example_alleles
    String cross
    File real_phases
    Float cmBymb
  }

  command <<<

        R --vanilla --no-save <<RSCRIPT
          library(supermassa4onemap)
          library(onemap)
          library(updog)
          library(reshape2)
          library(vcfR)
          library(doParallel)
          library(GUSMap)

          args = commandArgs(trailingOnly=TRUE)

          system("cp ~{freebayes_ref_depth} .")
          system("cp ~{freebayes_alt_depth} .")
          system("cp ~{gatk_ref_depth} .")
          system("cp ~{gatk_alt_depth} .")
          system("cp ~{gatk_example_alleles} .")
          system("cp ~{freebayes_example_alleles} .")
          method_name <- "~{methodName}"
          tot_mks_file <- "~{tot_mks}"
          simu_vcf_file <- "~{simu_vcf}"
	  max.cores <- 20
          vcf_file <- "~{vcf_file}"
          cross <- "~{cross}"
          cMbyMb <- ~{cmBymb}
          real_phases <- read.table("~{real_phases}")
          source("/opt/scripts/functions_simu.R")


         ## KNOWN VARIANTS
          tot_mks <- read.table(tot_mks_file)
          
          if(cross == "F1"){
            cross <- "outcross"
            f1 = NULL
          } else if (cross == "F2"){
            cross <- "f2 intercross"
            f1 = "F1"
          }
          
          # READING DATA FROM SIMULATED POPULATION
          simu <- read.vcfR(simu_vcf_file)
          gab <- onemap_read_vcfR(vcfR.object=simu,
                                  cross= cross,
                                  parent1="P1",
                                  parent2="P2",
                                  f1 = f1)
          
          ## READING FINAL VCF FROM PIPELINE
          vcf <- read.vcfR(vcf_file)
          df <- onemap_read_vcfR(vcfR.object=vcf,
                                 cross= cross,
                                 parent1="P1",
                                 parent2="P2",
                                 f1 = f1)
          
          ## FILTERS REPORT
          SNPcall <- method_name
          Genocall <- "df"
          CountsFrom <- "vcf"
          filters_tab <- create_filters_report(df, SNPcall, CountsFrom, Genocall)
          
          ## MAPS REPORT - DF
          
          ## Without false SNPs
          times <-system.time(create_maps_report(input.seq = filters_tab,
                                                 tot_mks = tot_mks, gab = gab,
                                                 SNPcall , Genocall,
                                                 fake= F, CountsFrom,cMbyMb))
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", FALSE)
          times <- data.frame(meth = outname, time = times[3])
          
          # With false SNPs
          times_temp <-system.time(create_maps_report(input.seq = filters_tab,
                                                      tot_mks = tot_mks, gab = gab,
                                                      SNPcall = method_name, Genocall= "df",
                                                      fake= T, CountsFrom="vcf",cMbyMb))
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", TRUE)
          times_temp <- data.frame(meth = outname, time = times_temp[3])
          times <- rbind(times, times_temp)
          
          errors_tab <- create_errors_report(onemap_obj = df, gab,
                                             SNPcall, Genocall,
                                             CountsFrom)
          
          # Maps report - df global error 0.05
          aval.df0.05 <- create_probs(onemap.obj = df, global_error = 0.05)
          
          Genocall <- "df0.05"
          filters_tab <- create_filters_report(aval.df0.05, SNPcall, CountsFrom, Genocall)
          
          fake <- F
          
          times_temp <- system.time(create_maps_report(filters_tab,
                                                       tot_mks, gab,
                                                       SNPcall, Genocall,
                                                       fake, CountsFrom,cMbyMb))
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
          times_temp <- data.frame(meth = outname, time = times_temp[3])
          times <- rbind(times, times_temp)
          
          fake <- T
          times_temp <- system.time(create_maps_report(filters_tab, tot_mks, gab,
                                                       SNPcall, Genocall,
                                                       fake, CountsFrom,cMbyMb))
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
          times_temp <- data.frame(meth = outname, time = times_temp[3])
          times <- rbind(times, times_temp)
          
          create_errors_report(aval.df0.05, gab,
                               SNPcall, Genocall,
                               CountsFrom)
          
          # MAPS REPORT - GQ
          aval.gq <- extract_depth(vcfR.object=vcf,
                                   onemap.object=df,
                                   vcf.par="GQ",
                                   parent1="P1",
                                   parent2="P2",
                                   f1 = f1,
                                   recovering=FALSE)
          
          aval.gq <- create_probs(onemap.obj = df, genotypes_errors=aval.gq)
          
          Genocall <- "GQ"
          filters_tab <- create_filters_report(aval.gq, SNPcall, CountsFrom, Genocall)
          
          fake <- F
          
          times_temp <- system.time(create_maps_report(filters_tab,
                                                       tot_mks, gab,
                                                       SNPcall, Genocall,
                                                       fake, CountsFrom,cMbyMb))
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
          times_temp <- data.frame(meth = outname, time = times_temp[3])
          times <- rbind(times, times_temp)
          
          fake <- T
          times_temp <- system.time(create_maps_report(filters_tab, tot_mks, gab,
                                                       SNPcall, Genocall,
                                                       fake, CountsFrom,cMbyMb))
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
          times_temp <- data.frame(meth = outname, time = times_temp[3])
          times <- rbind(times, times_temp)
          
          
          create_errors_report(aval.gq, gab,
                               SNPcall, Genocall,
                               CountsFrom)
          
          # OTHER TOOLS
          ## With depths from vcf
          
          updog.aval <- updog_genotype(
            vcfR.object=vcf,
            onemap.object=df,
            vcf.par="AD",
            parent1="P1",
            parent2="P2",
            f1 = f1,
            recovering=TRUE,
            mean_phred=20,
            cores=3,
            depths=NULL,
            global_error = NULL,
            use_genotypes_errors = FALSE,
            use_genotypes_probs = TRUE)
          
          updog0.05.aval <- create_probs(onemap.obj = updog.aval, global_error = 0.05)
          
          supermassa.aval <- supermassa4onemap::supermassa_genotype(
            vcfR.object=vcf,
            onemap.object = df,
            vcf.par = "AD",
            parent1 = "P1",
            parent2 = "P2",
            f1 = f1,
            recovering = TRUE,
            mean_phred = 20,
            cores = 3,
            depths = NULL,
            global_error = NULL,
            use_genotypes_errors = FALSE,
            use_genotypes_probs = TRUE)
          
          supermassa0.05.aval <- create_probs(onemap.obj = supermassa.aval, global_error = 0.05)
          
          polyrad.aval <- polyRAD_genotype(
            vcf=vcf_file,
            onemap.obj=df,
            parent1="P1",
            parent2="P2",
            f1 = f1,
            crosstype=cross,
            global_error = NULL,
            use_genotypes_errors = FALSE,
            use_genotypes_probs = TRUE)
          
          polyrad0.05.aval <- create_probs(onemap.obj = polyrad.aval, global_error = 0.05)
          
          metodologies <- list(updog = updog.aval, supermassa = supermassa.aval, polyrad = polyrad.aval,
                               updog0.05 = updog0.05.aval, supermassa0.05 = supermassa0.05.aval, polyrad0.05 = polyrad0.05.aval)
          for (metodology in names(metodologies)){
            error_aval <- metodologies[[metodology]]
            ## Filters
            Genocall <- metodology
            SNPcall <- method_name
            CountsFrom <- "vcf"
            
            filters_tab <- create_filters_report(error_aval, SNPcall, CountsFrom, Genocall)
            
            ## Maps
            fake <- F
            times_temp <- system.time(create_maps_report(input.seq = filters_tab,
                                                         tot_mks, gab,
                                                         SNPcall, Genocall,
                                                         fake, CountsFrom,cMbyMb))
            
            outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
            times_temp <- data.frame(meth = outname, time = times_temp[3])
            times <- rbind(times, times_temp)
            
            fake <- T
            times_temp <- system.time(create_maps_report(input.seq = filters_tab,
                                                         tot_mks, gab,
                                                         SNPcall, Genocall,
                                                         fake, CountsFrom,cMbyMb))
            
            outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
            times_temp <- data.frame(meth = outname, time = times_temp[3])
            times <- rbind(times, times_temp)
            
            ## Errors info
            
            create_errors_report(error_aval, gab,
                                 SNPcall, Genocall,
                                 CountsFrom)
            
          }
          
          ## Depths from bam
          depths.alt <- read.table(paste0(method_name, "_alt_depth_bam.txt"), header = T)
          depths.ref <- read.table(paste0(method_name, "_ref_depth_bam.txt"), header = T)
          
          depths <- list("ref"=depths.ref, "alt"=depths.alt)
          CountsFrom <- "bam"
          updog.aval.bam <- updog_genotype(
            vcfR.object=vcf,
            onemap.object=df,
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
          
          updog0.05.aval.bam <- create_probs(onemap.obj = updog.aval.bam, global_error = 0.05)
          
          supermassa.aval.bam <- supermassa_genotype(
            vcfR.object=vcf,
            onemap.object = df,
            vcf.par = "AD",
            parent1 = "P1",
            parent2 = "P2",
            f1 = f1,
            recovering = TRUE,
            mean_phred = 20,
            cores = 3,
            depths = depths,
            global_error = NULL,
            use_genotypes_errors = FALSE,
            use_genotypes_probs = TRUE)
          
          supermassa0.05.aval.bam <- create_probs(onemap.obj = supermassa.aval.bam, global_error = 0.05)
          
          if(tail(strsplit(vcf_file, "[.]")[[1]],1) =="gz") {
            vcf.temp <- paste0(method_name,".", sample(1000,1), ".vcf")
            system(paste0("zcat ", vcf_file, " > ", vcf.temp))
            vcf_file <- vcf.temp
          }
          
          new.vcf <- make_vcf(vcf_file, depths, method_name)
          
          polyrad.aval.bam <- polyRAD_genotype(
            vcf=new.vcf,
            onemap.obj=df,
            parent1="P1",
            parent2="P2",
            f1 = f1,
            crosstype=cross,
            global_error = NULL,
            use_genotypes_errors = FALSE,
            use_genotypes_probs = TRUE)
          
          polyrad0.05.aval.bam <- create_probs(onemap.obj = polyrad.aval.bam, global_error = 0.05)
          
          metodologies <- list(updog = updog.aval.bam, supermassa= supermassa.aval.bam, polyrad=polyrad.aval.bam,
                               updog0.05 = updog0.05.aval.bam, supermassa0.05= supermassa0.05.aval.bam, polyrad0.05=polyrad0.05.aval.bam)
          for (metodology in names(metodologies)){
            error_aval <- metodologies[[metodology]]
            ## Filters
            Genocall <- metodology
            CountsFrom <- "bam"
            
            filters_tab <- create_filters_report(error_aval, SNPcall, CountsFrom, Genocall)
            
            ## Maps
            fake <- F
            times_temp <- system.time(create_maps_report(filters_tab,
                                                         tot_mks, gab,
                                                         SNPcall, Genocall,
                                                         fake, CountsFrom,cMbyMb))
            
            outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
            times_temp <- data.frame(meth = outname, time = times_temp[3])
            times <- rbind(times, times_temp)
            
            fake = T
            times_temp <- system.time(create_maps_report(filters_tab,
                                                         tot_mks, gab,
                                                         SNPcall, Genocall,
                                                         fake, CountsFrom,cMbyMb))
            
            outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
            times_temp <- data.frame(meth = outname, time = times_temp[3])
            times <- rbind(times, times_temp)
            
            ## Errors info
            errors_tab <- create_errors_report(onemap_obj = error_aval, gab = gab,
                                               SNPcall, Genocall,
                                               CountsFrom)
          }
          
          ## Gusmap maps
          Genocall <- "gusmap"
          
          fake <- F
          CountsFrom <- "vcf"
          times_temp <- system.time(create_gusmap_report(vcf_file, gab,SNPcall, Genocall,
                                                         fake, CountsFrom))
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
          times_temp <- data.frame(meth = outname, time = times_temp[3])
          times <- rbind(times, times_temp)
          
          fake <- T
          times_temp <- system.time(create_gusmap_report(vcf_file, gab,SNPcall, Genocall,
                                                         fake, CountsFrom))
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
          times_temp <- data.frame(meth = outname, time = times_temp[3])
          times <- rbind(times, times_temp)
          
          CountsFrom <- "bam"
          fake <- F
          times_temp <- system.time(create_gusmap_report(new.vcf, gab,SNPcall, Genocall,
                                                         fake, CountsFrom))
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
          times_temp <- data.frame(meth = outname, time = times_temp[3])
          times <- rbind(times, times_temp)
          
          fake <- T
          times_temp <- system.time(create_gusmap_report(new.vcf, gab,SNPcall, Genocall,
                                                         fake, CountsFrom))
          
          outname <- paste0("map_", SNPcall, "_", CountsFrom, "_", Genocall, "_", fake)
          times_temp <- data.frame(meth = outname, time = times_temp[3])
          times <- rbind(times, times_temp)
          
          write.table(times, paste0(method_name,"_times.txt"))
          
          Genocall <- c("df", "GQ", "updog", "supermassa", "polyrad", "gusmap",
                        "df0.05", "updog0.05", "supermassa0.05", "polyrad0.05")
          fake <- c(TRUE, FALSE)
          CountsFrom <- c("vcf", "bam")
          
          all_maps <- all_errors <- all_filters <- data.frame()
          all_RDatas <- list()
          z <- 1
          names_RDatas <- vector()
          for(i in 1:length(Genocall)){
            for(j in 1:length(CountsFrom)){
              for(w in 1:length(fake)){
                if(CountsFrom[j] == "bam" & (Genocall[i] == "df" | Genocall[i] == "GQ" |  Genocall[i] == "df0.05" )){
                } else {
                  cat(paste0("map_",method_name,"_",CountsFrom[j], "_", Genocall[i], "_",fake[w], ".txt"), "\n")
                  temp_map <- read.table(paste0("map_",method_name,"_",CountsFrom[j], "_", Genocall[i], "_",fake[w], ".txt"), header=T)
                  load(paste0("map_",method_name,"_",CountsFrom[j], "_", Genocall[i], "_",fake[w], ".RData"))
                  all_RDatas[[z]] <- map_df
                  names_RDatas <- c(names_RDatas, paste0("map_",method_name,"_",CountsFrom[j], "_", Genocall[i], "_",fake[w]))
                  z <- z+1
                  all_maps <- rbind(all_maps, temp_map)
                  if(Genocall[i] == "gusmap"){
                    
                  } else {
                    cat(paste0("errors_",method_name,"_",CountsFrom[j], "_", Genocall[i], ".txt"), "\n")
                    temp_error <- read.table(paste0("errors_",method_name,"_",CountsFrom[j], "_", Genocall[i],".txt"), header=T)
                    all_errors <- rbind(all_errors, temp_error)
                    cat(paste0("filters_",method_name,"_",CountsFrom[j], "_", Genocall[i],".txt"), "\n")
                    temp_filters <- read.table(paste0("filters_",method_name,"_",CountsFrom[j], "_", Genocall[i],".txt"), header=T)
                    all_filters <- rbind(all_filters, temp_filters)
                  }
                }
              }
            }
          }
          
          names(all_RDatas) <- names_RDatas
          write.table(all_maps, paste0(method_name,"_maps.txt"))
          write.table(all_errors, paste0(method_name,"_errors.txt"))
          write.table(all_filters, paste0(method_name,"_filters.txt"))
          save(all_RDatas, file= paste0(method_name, "_RDatas.RData"))

        RSCRIPT
  >>>

  runtime{
    docker:"taniguti/onemap"
    mem:"--nodes=1"
    cpu:20
    time:"192:00:00"
  }

  output{
    File times = "~{methodName}_times.txt"
    File all_maps = "~{methodName}_maps.txt"
    File all_errors = "~{methodName}_errors.txt"
    File all_filters = "~{methodName}_filters.txt"
    File all_RDatas = "~{methodName}_RDatas.RData"
  }
}

task CreateTables{
  input{
    Int depth
    Int seed
    File tot_mks
    File gatk_ref_depth
    File gatk_ref_depth_bam
    File gatk_alt_depth
    File gatk_alt_depth_bam
    File freebayes_ref_depth_bam
    File freebayes_alt_depth_bam
    File freebayes_ref_depth
    File freebayes_alt_depth
    Array[File] all_maps
    Array[File] all_errors
    Array[File] all_filters
    Array[File] times
    Array[File] all_RDatas
  }

  command <<<

        R --vanilla --no-save <<RSCRIPT

          library(reshape2)
          library(vcfR)
          system("cp ~{sep= " " times } .")
          system("cp ~{sep= " " all_maps } .")
          system("cp ~{sep= " " all_errors } .")
          system("cp ~{sep= " " all_filters } .")
          system("cp ~{sep= " " all_RDatas } .")

          tot_mks <- read.table("~{tot_mks}")
          method <- c("gatk", "freebayes")
          meth.bam <-  c("polyrad", "updog", "supermassa")
          meth.filt <- c("dfAndGQ", "polyrad", "updog", "supermassa")
          meth.geno <- c("df", "GQ", "polyrad", "updog", "supermassa")
          seed <- ~{seed}
          depth <- ~{depth}

          # Functions

           joint_depths <- function(depth_matrix_alt, depth_matrix_ref, CountsFrom, SNPcall, depth,seed, genotype_meth){
              depth_matrix <- list(depth_matrix_alt, depth_matrix_ref)
              alle_name <- c("alt", "ref")
              allele <- list()
              for(i in 1:length(depth_matrix)){
                depth2 <- read.table(depth_matrix[[i]], header = T, stringsAsFactors = F)
                depth3 <- as.data.frame(apply(depth2, 2, as.integer))
                depth3 <- cbind(MKS= rownames(depth2), depth3)
                depth3 <- depth3[,sort(colnames(depth3))]
                allele[[i]] <- melt(depth3)
                colnames(allele[[i]]) <- c("mks", "ind", paste0(alle_name[i]))
              }
              alleles <- merge(allele[[1]], allele[[2]])
              alleles <- cbind(seed= seed, depth = depth, "SNPcall" = SNPcall, "CountsFrom" = CountsFrom, alleles)
              alleles[,4] <- as.character(alleles[,4])
              alleles[,5] <- as.character(alleles[,5])
              return(alleles)
            }


            ########################################################################################
            # Table1: GenoCall; mks; ind; SNPcall; CountsFrom; alt; ref; gabGT; methGT; A; AB; BA; B
            ########################################################################################
            df_tot_all <- times <- maps_tot <- filters_tot <- vector()
            tot_RDatas <- list()
            for(i in 1:length(method)){

              if(method[i] == "gatk"){
                alt.depth.bam <- "~{gatk_alt_depth_bam}"
                ref.depth.bam <- "~{gatk_ref_depth_bam}"
                alt.depth <- "~{gatk_alt_depth}"
                ref.depth <- "~{gatk_ref_depth}"
              } else{
                alt.depth.bam <- "~{freebayes_alt_depth_bam}"
                ref.depth.bam <- "~{freebayes_ref_depth_bam}"
                alt.depth <- "~{freebayes_alt_depth}"
                ref.depth <- "~{freebayes_ref_depth}"
              }

              ## Depths by bam

              df_tot_bam <- joint_depths(depth_matrix_alt = alt.depth.bam, depth_matrix_ref = ref.depth.bam,
                                         CountsFrom = "bam", SNPcall = method[i], depth = depth, seed = seed, genotype_meth = meth.bam)

              ## Depths by softwares
              df_tot <- joint_depths(alt.depth,ref.depth, "vcf", method[i], depth, seed, meth.geno)

              df_tot_tot <- rbind(df_tot_bam, df_tot)

              chr <- unique(sapply(strsplit(df_tot_tot[,5], "_"), "[",1))
              all_errors <- read.table(paste0(method[i],"_errors.txt"), header = T)
              all_errors[,5] <- paste0(chr, "_",all_errors[,5])
              colnames(all_errors) <- c("SNPcall", "Genocall", "CountsFrom", "ind", "mks", "gabGT", "methGT", "A", "AB", "BA", "B")
              df_meth <- merge(all_errors, df_tot_tot, by = c("SNPcall", "CountsFrom", "ind", "mks"))

              df_tot_all <- rbind(df_tot_all, df_meth)

              ########################################################
              # Table2: seed; CountsFrom; ErrorProb; SNPcall; MK; rf; phases; real_phases
              ########################################################

              map_df <- read.table(paste0(method[i], "_maps.txt"), header = T)
              map_df <- data.frame("seed" = seed, "depth" = depth, map_df)

              maps_tot <- rbind(maps_tot, map_df)

              ##########################################################################
              # Table3: CountsFrom; seed; SNPcall; GenoCall; n_mks; distorted; redundant
              ##########################################################################

              filters_temp <- read.table(paste0(method[i], "_filters.txt"), header = T)
              filters_temp <- data.frame("seed" = seed, "depth" = depth, filters_temp)

              filters_tot <- rbind(filters_tot, filters_temp)

              ###########################################################################
              # Table4: CountsFrom; seed; SNPcall; GenoCall
              ###########################################################################

              CountsFrom <- vector()
              times_temp <- read.table(paste0(method[i], "_times.txt"), stringsAsFactors = F)
              temp <- strsplit(times_temp[,1], "_")
              temp <- do.call(rbind, temp)
              CountsFrom[which(temp[,3] == "bam")] <- "bam"
              CountsFrom[which(temp[,3] != "bam")] <- "vcf"
              SNPcall <- temp[,2]
              temp[which(temp[,3] == "bam"),3] <- temp[which(temp[,3] == "bam"),4]
              Genocall <- temp[,4]
              real.mks <- temp[,5]
              times <- rbind(times, data.frame(CountsFrom, seed, depth, SNPcall, Genocall, real.mks, times = times_temp[,2]))

              ###########################################################################
              # Table6: list of RDatas with name CountsFrom; seed; SNPcall; GenoCall
              ###########################################################################

              load(paste0(method[i], "_RDatas.RData"))
              tot_RDatas <- c(tot_RDatas, all_RDatas)

            }

            names_RDatas <- names(tot_RDatas)
            new_names <- paste0(seed, "_", depth, "_", names_RDatas)
            names(tot_RDatas) <- new_names

            saveRDS(df_tot_all, file = "data1_depths_geno_prob.rds")
            saveRDS(maps_tot, file = "data2_maps.rds")
            saveRDS(filters_tot, file = "data3_filters.rds")
            saveRDS(times, file= "data4_times.rds")
            save(tot_RDatas, file = "data6_RDatas.RData")

      RSCRIPT
  >>>

  runtime{
    docker:"taniguti/onemap"
    mem:"--nodes=1"
    cpu:1
    time:"01:00:00"
  }

  output{
    File data1_depths_geno_prob = "data1_depths_geno_prob.rds"
    File data2_maps = "data2_maps.rds"
    File data3_filters = "data3_filters.rds"
    File data4_times   = "data4_times.rds"
    File data6_RDatas  = "data6_RDatas.RData"
  }
}

