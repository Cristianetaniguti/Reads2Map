version 1.0

task vcf2onemap {
   input {
     File vcf_file
     String cross
     String parent1
     String parent2
   }

   Int disk_size = ceil(size(vcf_file, "GiB") * 2)
   Int memory_size = ceil(size(vcf_file, "MiB") * 2 + 2000)

   command <<<

        R --vanilla --no-save <<RSCRIPT
          library(onemap)
          library(vcfR)

          cross <- "~{cross}"

          if(cross == "F1"){
            cross <- "outcross"
            f1 = NULL
          } else if (cross == "F2"){
            cross <- "f2 intercross"
            f1 = "F1"
          }

          ## READING VCF FROM PIPELINE
          vcfR.obj <- read.vcfR("~{vcf_file}")
          save(vcfR.obj, file="vcfR_obj.RData")

          onemap.obj <- onemap_read_vcfR(vcfR.object= vcfR.obj,
                                         cross= cross,
                                         parent1="~{parent1}",
                                         parent2="~{parent2}",
                                         f1 = f1)

          save(onemap.obj, file=paste0("vcf_onemap.obj.RData"))

        RSCRIPT

    >>>

    runtime {
      docker:"cristaniguti/reads2map:0.0.4"
      cpu:1
      # Cloud
      memory:"~{memory_size} MiB"
      disks:"local-disk " + disk_size + " HDD"
      # Slurm
      job_name: "vcf2onemap"
      mem:"~{memory_size}M"
      time:"10:00:00"
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Convert VCF file to onemap object. See [OneMap](https://github.com/Cristianetaniguti/onemap) for more information."
    }

    output {
      File onemap_obj = "vcf_onemap.obj.RData"
      File vcfR_obj = "vcfR_obj.RData"
    }
}

task FiltersReport {
  input {
    File onemap_obj
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    Int seed
    Int depth
  }

  Int disk_size = ceil(size(onemap_obj, "GiB") * 2)
  Int memory_size = ceil(size(onemap_obj, "MiB") * 3 + 2000)

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)
      library(Reads2MapTools)

      temp <- load("~{onemap_obj}")
      temp.obj <- get(temp)
      onemap_obj_filtered <- create_filters_report_simu(temp.obj, "~{SNPCall_program}",
                                                  "~{CountsFrom}", "~{GenotypeCall_program}",
                                                   ~{seed}, ~{depth}, threshold = NULL) # Threshold define the genotype probability filter

      save(onemap_obj_filtered, file="onemap_obj_filtered.RData")

    RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/reads2map:0.0.4"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "Filters"
    mem:"~{memory_size}M"
    time:"05:00:00"
  }

  meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Filter simulated markers in onemap object by missing data, segregation distortion and redundant markers. See [Reads2MapTools](https://github.com/Cristianetaniguti/Reads2MapTools) for more information."
  }

  output {
    File filters_report = "~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}_~{seed}_~{depth}_filters_report.tsv.gz"
    File onemap_obj_filtered = "onemap_obj_filtered.RData"
  }
}

task FiltersReportEmp {
  input {
    File onemap_obj
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    String chromosome
  }

  Int disk_size = ceil(size(onemap_obj, "GiB") * 2)
  Int memory_size = ceil(size(onemap_obj, "MiB") * 2 + 3000)

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)
      library(Reads2MapTools)

      temp <- load("~{onemap_obj}")
      temp.obj <- get(temp)
      onemap_obj_filtered <- create_filters_report_emp(temp.obj, "~{SNPCall_program}",
                                           "~{CountsFrom}", "~{GenotypeCall_program}",
                                           "~{chromosome}", threshold = NULL)
      save(onemap_obj_filtered, file="onemap_obj_filtered.RData")

    RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/reads2map:0.0.4"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "Filters"
    mem:"~{memory_size}M"
    time:"05:00:00"
  }

  meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Filter empirical markers in onemap object by missing data, segregation distortion and redundant markers. See [Reads2MapTools](https://github.com/Cristianetaniguti/Reads2MapTools) for more information."
  }

  output {
    File filters_report = "~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}_filters_report.tsv.gz"
    File onemap_obj_filtered = "onemap_obj_filtered.RData"
  }
}

task MapsReport {
  input {
   File onemap_obj
   File ref_alt_alleles
   File simu_onemap_obj
   String SNPCall_program
   String GenotypeCall_program
   String CountsFrom
   File simulated_phases
   Int seed
   Int depth
   Int max_cores
  }

  Int disk_size = ceil((size(onemap_obj, "GiB") * 2) + size(ref_alt_alleles, "GiB") + (size(simu_onemap_obj, "GiB") * 2) + size(simulated_phases, "GiB"))
  Int memory_size = ceil(size(onemap_obj, "MiB") * 3 + 4000)

  command <<<
      R --vanilla --no-save <<RSCRIPT
      library(onemap)
      library(Reads2MapTools)

      if(~{max_cores} > 4) cores = 4 else cores = ~{max_cores}

      filtered_onemap <- load("~{onemap_obj}")
      filtered_onemap <- get(filtered_onemap)

      simu_onemap_obj <- load("~{simu_onemap_obj}")
      simu_onemap_obj <- get(simu_onemap_obj)
      ref_alt_alleles <- read.table("~{ref_alt_alleles}")
      simulated_phases <- read.table("~{simulated_phases}")

      ## With false SNPs
      times_fake <-system.time(info_fake <- create_maps_report_simu(input.seq = filtered_onemap,
                                                  tot_mks = ref_alt_alleles, gab = simu_onemap_obj,
                                                  "~{SNPCall_program}" , "~{GenotypeCall_program}",
                                                  fake= "with-false", "~{CountsFrom}", simulated_phases,
                                                  ~{seed}, ~{depth}, cores))

      times <- data.frame(seed = ~{seed}, depth = ~{depth}, SNPCall = "~{SNPCall_program}",
                          CountsFrom = "~{CountsFrom}", GenoCall =  "~{GenotypeCall_program}", fake = "with-false",
                          time = times_fake[3])

      # It will not run if all markers are true markers
      if(all(info_fake[[2]][,"real.mks"] == "true marker")){
        cat("skip :) \n")
        times_temp <- times_fake
        info_correct <- update_fake_info(info_fake, simu_onemap_obj, ref_alt_alleles, simulated_phases)

      } else {
        ## Without false SNPs
        times_temp <-system.time(info_correct <- create_maps_report_simu(input.seq = filtered_onemap,
                                              tot_mks = ref_alt_alleles, gab = simu_onemap_obj,
                                              "~{SNPCall_program}" , "~{GenotypeCall_program}",
                                              fake= "without-false", "~{CountsFrom}", simulated_phases,
                                              ~{seed}, ~{depth}, cores))
      }

      # Joint maps data.frames
      map_joint <- rbind(info_fake[[2]], info_correct[[2]])
      vroom::vroom_write(map_joint, "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_~{seed}_~{depth}_map_report.tsv.gz", num_threads = ~{max_cores})

      # Joint RDatas
      RDatas_joint <- list()
      RDatas_joint[[1]] <- info_fake[[1]]
      RDatas_joint[[2]] <- info_correct[[1]]
      names(RDatas_joint) <- c("map_~{seed}_~{depth}_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_TRUE",
                               "map_~{seed}_~{depth}_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_FALSE")
      save(RDatas_joint, file= "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_~{seed}_~{depth}.RData")

      # Joint times data.frames
      times_temp <- data.frame(seed = ~{seed}, depth = ~{depth}, SNPCall = "~{SNPCall_program}",
                               CountsFrom = "~{CountsFrom}", GenoCall =  "~{GenotypeCall_program}", fake = "without-false",
                               time = times_temp[3])

      times <- rbind(times, times_temp)
      vroom::vroom_write(times,"~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_~{seed}_~{depth}_times_report.tsv.gz", num_threads = ~{max_cores})

      RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/reads2map:0.0.4"
    cpu:4
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "MapsReport"
    mem:"~{memory_size}M"
    time:"48:00:00"
  }

  meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Estimate genetic distances by OneMap HMM multi-point approach in a set o simulated markers ordered by genomic position with and without false-positives. See [Reads2MapTools](https://github.com/Cristianetaniguti/Reads2MapTools) for more information."
  }

  output {
    File maps_report = "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_~{seed}_~{depth}_map_report.tsv.gz"
    File maps_RData = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_~{seed}_~{depth}.RData"
    File times = "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_~{seed}_~{depth}_times_report.tsv.gz"
  }
}

task ErrorsReport {
  input {
    File onemap_obj
    File vcfR_obj
    File simu_vcfR
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    Int max_cores
    Int seed
    Int depth
  }

  Int disk_size = ceil((size(onemap_obj, "GiB") * 2) + size(vcfR_obj, "GiB") + (size(simu_vcfR, "GiB") * 2))
  Int memory_size = ceil(size(onemap_obj, "MiB") * 3 + 3000)

  command <<<
      R --vanilla --no-save <<RSCRIPT

        library(onemap)
        library(tidyverse)

        temp <- load("~{onemap_obj}")
        df <- get(temp)

        temp <- load("~{vcfR_obj}")
        vcf <- get(temp)

        p <- create_depths_profile(onemap.obj = df, vcfR.object = vcf, parent1 = "P1",
                                   parent2 = "P2", vcf.par = "AD",recovering = FALSE,
                                   GTfrom = "onemap", alpha=0.1,
                                   rds.file = paste0("vcf_depths.rds"))

        df <- readRDS(paste0("vcf_depths.rds"))
        df <- cbind(seed = "~{seed}"  , depth = "~{depth}", SNPCall = "~{SNPCall_program}", CountsFrom = "~{CountsFrom}",
                    GenoCall="~{GenotypeCall_program}", df)

        simu <- load("~{simu_vcfR}")
        vcf_simu <- get(simu)

        gt.simu <- vcf_simu@gt[,-1]
        gt.simu <- as.data.frame(cbind(mks = vcf_simu@fix[,3], gt.simu))

        gt.simu <- data.frame(lapply(gt.simu, function(x) {
            x %>% gsub("[|]", "/", .) %>%
            gsub("[.]","missing", .) %>%
            gsub("0/0","homozygous-ref", .) %>%
            gsub("1/1","homozygous-alt", .) %>%
            gsub("0/1","heterozygous", .) %>%
            gsub("1/0","heterozygous", .)
        }))

        dptot <- gt.simu %>%
          pivot_longer(!mks, names_to = "ind", values_to = "gabGT") %>% inner_join(df)

        idx <- match(c("A", "AB", "BA", "B"), colnames(dptot))
        dptot <- cbind(dptot, errors = apply(dptot[,idx], 1, function(x) 1 - max(x)))

        vroom::vroom_write(dptot, "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_~{seed}_~{depth}_errors_report.tsv.gz", num_threads = ~{max_cores})

      RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/reads2map:0.0.4"
    cpu: max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "ErrorsReport"
    mem:"~{memory_size}M"
    time:"05:00:00"
  }

  meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Creates data.frame with probabilities used in the OneMap HMM. See [OneMap](https://github.com/Cristianetaniguti/onemap) for more information."
  }

  output {
    File errors_report = "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_~{seed}_~{depth}_errors_report.tsv.gz"
  }
}

task CheckDepths {
  input {
    File onemap_obj
    File vcfR_obj
    String parent1
    String parent2
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    Int? seed
    Int? depth
    Int max_cores
  }

  Int disk_size = ceil((size(onemap_obj, "GiB") * 2) + size(vcfR_obj, "GiB"))
  Int memory_size = ceil(size(onemap_obj, "MiB") * 3 + 3000)

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)

      temp <- load("~{onemap_obj}")
      df <- get(temp)

      temp <- load("~{vcfR_obj}")
      vcf <- get(temp)

      p <- create_depths_profile(onemap.obj = df, vcfR.object = vcf, parent1 = "~{parent1}",
      parent2 = "~{parent2}", vcf.par = "AD",recovering = FALSE, GTfrom = "vcf", alpha=0.1,
      rds.file = paste0("vcf_depths.rds"))

      df <- readRDS(paste0("vcf_depths.rds"))
      df <- cbind(SNPCall = "~{SNPCall_program}", CountsFrom = "~{CountsFrom}",
                  GenoCall="~{GenotypeCall_program}", df ~{", seed=" + seed} ~{", depth= " + depth})

      vroom::vroom_write(df, "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_errors_report.tsv.gz", num_threads = ~{max_cores})

    RSCRIPT

  >>>

  runtime {
    docker:"cristaniguti/reads2map:0.0.4"
    cpu: max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "CheckDepths"
    mem:"~{memory_size}M"
    time:"10:00:00"
  }

  meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Creates data.frame with reads depth information. See [OneMap](https://github.com/Cristianetaniguti/onemap) for more information."
  }

  output {
    File errors_report = "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_errors_report.tsv.gz"
  }

}

task MapsReportEmp {
  input {
   File sequence_obj
   String SNPCall_program
   String GenotypeCall_program
   String CountsFrom
   Int max_cores
  }

  Int disk_size = ceil((size(sequence_obj, "GiB") * 2))
  Int memory_size = ceil(size(sequence_obj, "MiB") * 3 + 4000)

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)
      library(Reads2MapTools)

      temp <- load("~{sequence_obj}")
      sequence <- get(temp)

      if(~{max_cores} > 4) cores = 4 else cores = ~{max_cores}

      times_temp <- system.time(df <- create_map_report_emp(input.seq = sequence, CountsFrom = "~{CountsFrom}",
                                      SNPCall = "~{SNPCall_program}", GenoCall="~{GenotypeCall_program}", max_cores = cores))

      vroom::vroom_write(df[[2]], "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_map_report.tsv.gz", num_threads = ~{max_cores})
      map_out <- df[[1]]
      save(map_out,  file = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.RData")

      times <- data.frame(SNPCall = "~{SNPCall_program}",
                          CountsFrom = "~{CountsFrom}",
                          GenoCall =  "~{GenotypeCall_program}",
                          time = times_temp[3])

      vroom::vroom_write(times, "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_times_report.tsv.gz", num_threads = ~{max_cores})

    RSCRIPT
  >>>

  runtime {
    docker:"cristaniguti/reads2map:0.0.4"
    cpu:max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "MapsEmp"
    mem:"~{memory_size}M"
    time:"48:00:00"
  }

  meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Estimate genetic distances by OneMap HMM multi-point approach in a set o empirical markers ordered by genomic position. See [Reads2MapTools](https://github.com/Cristianetaniguti/Reads2MapTools) for more information."
  }

  output {
    File maps_report = "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_map_report.tsv.gz"
    File maps_RData = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.RData"
    File times = "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_times_report.tsv.gz"
  }
}

# Genotype calling with updog, polyrad and supermassa
# Exclusive for biallelic markers, input VCF file are normalized with -m-any
task ReGenotyping {
  input {
    String GenotypeCall_program
    File vcf_file
    String cross
    String parent1
    String parent2
    Int max_cores
    Int ploidy
  }

  Int disk_size = ceil((size(vcf_file, "GiB") * 4))
  Int memory_size = ceil(size(vcf_file, "MiB") * 3 + 4000)

  command <<<
     R --vanilla --no-save <<RSCRIPT
       library(onemap)
       library(Reads2MapTools)

       method <- "~{GenotypeCall_program}"
       cross <- "~{cross}"

       if(cross == "F1"){
        cross <- "outcross"
       } else {
        stop("Invalid crosstype")
       }

       out_vcf <- "regeno.vcf.gz"

        if (method == "updog") {
          updog_genotype_vcf(vcf="~{vcf_file}",
                             vcf.par="AD",
                             out_vcf = out_vcf,
                             parent1="~{parent1}",
                             parent2="~{parent2}",
                             crosstype= cross,
                             cores=~{max_cores},
                             ploidy = ~{ploidy})
        } else if (method == "supermassa") {
          supermassa_genotype_vcf(vcf="~{vcf_file}",
                                  vcf.par="AD",
                                  out_vcf = out_vcf,
                                  parent1="~{parent1}",
                                  parent2="~{parent2}",
                                  crosstype= cross,
                                  cores=~{max_cores},
                                  ploidy = ~{ploidy})
        } else if (method == "polyrad") {
          ex <- strsplit(basename("~{vcf_file}"), split="[.]")[[1]]
          if(ex[length(ex)] == "gz") {
              system("gunzip -f ~{vcf_file}")
              vcf_in <- paste0(dirname("~{vcf_file}"), "/", paste0(ex[-length(ex)], collapse = "."))
          } else vcf_in <- "~{vcf_file}"

          polyRAD_genotype_vcf(vcf=vcf_in,
                               parent1="~{parent1}",
                               parent2="~{parent2}",
                               ploidy = ~{ploidy},
                               out_vcf = out_vcf)                                                 
        }

        system("gunzip regeno.vcf.gz")


     RSCRIPT
  >>>

  runtime {
    docker:"cristaniguti/reads2map:0.0.5"
    cpu: max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "ReGenotyping"
    mem:"~{memory_size}M"
    time:"10:00:00"
  }

  meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Run updog, polyRAD and SuperMASSA genotyping for F1 population. See [Reads2MapTools](https://github.com/Cristianetaniguti/Reads2MapTools) for more information."
  }

  output {
    File regeno_vcf = "regeno.vcf"
  }
}

task SetProbs {
  input {
    File vcf_file
    String cross
    String parent1
    String parent2
    String multiallelics
    String SNPCall_program
  }

  Int disk_size = ceil(size(vcf_file, "GiB") * 3)
  Int memory_size = ceil(size(vcf_file, "MiB") * 2 + 4000)

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)
      library(vcfR)

      cross <- "~{cross}"

      if(cross == "F1"){
         cross <- "outcross"
         f1 = NULL
      } else if (cross == "F2"){
         f1 = "F1"
      }

      vcf <- read.vcfR("~{vcf_file}")
      save(vcf, file = "vcfR.RData")

      if("~{multiallelics}") only_biallelic = FALSE else only_biallelic = TRUE

      onemap.obj <- onemap_read_vcfR(vcfR.object = vcf,
                                     cross= cross,
                                     parent1="~{parent1}",
                                     parent2="~{parent2}",
                                     f1 = f1, only_biallelic = only_biallelic)

      # if("~{SNPCall_program}" == "freebayes") par <- "GL" else par <- "PL"

      probs <- extract_depth(vcfR.object=vcf,
                               onemap.object=onemap.obj,
                               vcf.par= "GQ",
                               parent1="~{parent1}",
                               parent2="~{parent2}",
                               f1 = f1,
                               recovering=FALSE)

      probs_onemap_obj <- create_probs(input.obj = onemap.obj, genotypes_errors=probs)

      # If filter by genotype probability
      # onemap_prob <- filter_prob(probs_onemap_obj, threshold = threshold)
      # onemap_mis <- filter_missing(onemap_prob, threshold = 0.25)
      # globalerror_onemap_obj <- create_probs(input.obj = onemap_mis, global_error = 0.05)

      globalerror_onemap_obj <- create_probs(input.obj = onemap.obj, global_error = 0.05)

      save(probs_onemap_obj, file="probs_onemap_obj.RData")
      save(globalerror_onemap_obj, file="globalerror_onemap_obj.RData")

    RSCRIPT

  >>>
  runtime {
    docker:"cristaniguti/reads2map:0.0.4"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "SetProbs"
    mem:"~{memory_size}M"
    time:"10:00:00"
  }

  meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description:"Get GQ or PL values from VCF file and set them as genotype probabilities in OneMap object. Also produce an OneMap object with global error rate of 0.05. See [OneMap](https://github.com/Cristianetaniguti/onemap) for more information."
  }

  output {
    File probs_onemap_obj = "probs_onemap_obj.RData"
    File globalerror_onemap_obj = "globalerror_onemap_obj.RData"
    File vcfR_obj = "vcfR.RData"
  }
}


task SetProbsDefault {
  input {
    File vcf_file
    File? multiallelics_mchap
    String mchap
    String SNPCall_program
    String cross
    String parent1
    String parent2
    String multiallelics
  }

  Int disk_size = ceil(size(vcf_file, "GiB") * 3)
  Int memory_size = ceil(size(vcf_file, "MiB") * 2 + 3000)

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)
      library(vcfR)

      cross <- "~{cross}"

      if(cross == "F1"){
        cross <- "outcross"
        f1 = NULL
      } else if (cross == "F2"){
        f1 = "F1"
      }

      if(as.logical("~{mchap}") & "~{SNPCall_program}" == "gatk") vcf <- read.vcfR("~{multiallelics_mchap}") else  vcf <- read.vcfR("~{vcf_file}")
      save(vcf, file = "vcfR.RData")

      if("~{multiallelics}") only_biallelic = FALSE else only_biallelic = TRUE

      onemap.obj <- onemap_read_vcfR(vcfR.object = vcf,
                                     cross= cross,
                                     parent1="~{parent1}",
                                     parent2="~{parent2}",
                                     f1 = f1, only_biallelic = only_biallelic)

      # if("~{SNPCall_program}" == "freebayes") par <- "GL" else par <- "PL"

      probs <- extract_depth(vcfR.object=vcf,
                               onemap.object=onemap.obj,
                               vcf.par= "GQ",
                               parent1="~{parent1}",
                               parent2="~{parent2}",
                               f1 = f1,
                               recovering=FALSE)

      probs_onemap_obj <- create_probs(input.obj = onemap.obj, genotypes_errors=probs)

      # onemap_prob <- filter_prob(probs_onemap_obj, threshold = threshold)
      # onemap_mis <- filter_missing(onemap_prob, threshold = 0.25)
      # globalerror_onemap_obj <- create_probs(input.obj = onemap_mis, global_error = 0.05)
      globalerror_onemap_obj <- create_probs(input.obj = onemap.obj, global_error = 0.05)

      default_onemap_obj <- create_probs(input.obj = onemap.obj, global_error = 10^(-5))
      save(default_onemap_obj, file="default_onemap_obj.RData")
      save(probs_onemap_obj, file="probs_onemap_obj.RData")
      save(globalerror_onemap_obj, file="globalerror_onemap_obj.RData")

    RSCRIPT

  >>>
  runtime {
    docker:"cristaniguti/reads2map:0.0.4"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "SetProbsDefault"
    mem:"~{memory_size}M"
    time:"10:00:00"
  }

  meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Set OneMap object genotype probabilities to value used in OneMap < 3.0 and to 0.05. See [OneMap](https://github.com/Cristianetaniguti/onemap) for more information."
  }

  output {
    File probs_onemap_obj = "probs_onemap_obj.RData"
    File globalerror_onemap_obj = "globalerror_onemap_obj.RData"
    File default_onemap_obj = "default_onemap_obj.RData"
    File vcfR_obj = "vcfR.RData"
  }
}

task RemoveNonInformative {
  input {
    File vcf_file
    String parent1
    String parent2
    String replaceADbyMissing
  }

  Int disk_size = ceil(size(vcf_file, "GiB") * 2)
  Int memory_size = 3000

  command <<<
      R --vanilla --no-save <<RSCRIPT
        library(Reads2MapTools)
        # Replace AD and DP by 0 when GT is missing
        # Remove non-informatives
        remove_non_informative("~{vcf_file}", 
                                P1 = "~{parent1}", 
                                P2 = "~{parent2}",
                                replaceAD = "~{replaceADbyMissing}",
                                out.vcf = "filtered.vcf.gz")

      RSCRIPT
  >>>

  runtime {
      docker:"cristaniguti/reads2map:0.0.4"
      cpu:1
      # Cloud
      memory:"~{memory_size} MiB"
      disks:"local-disk " + disk_size + " HDD"
      # Slurm
      job_name: "FilterSegregation"
      mem:"~{memory_size}M"
      time:"10:00:00"
  }

  meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Filter VCF file according to markers segregation information and replace DP/AD VCF information by 0 when GT is missing. See [Reads2MapTools](https://github.com/Cristianetaniguti/Reads2MapTools) for more information."
  }

  output {
    File vcf_filtered = "filtered.vcf.gz"
  }

}

task QualPlots {
    input {
        File Total
    }

    Int disk_size = ceil(size(Total, "GB") + 1)
    Int memory_size = 3000

    command <<<
        R --vanilla --no-save <<RSCRIPT
            library(ggplot2)
            library(dplyr)
            library(tidyr)

            tot <- read.table("~{Total}", header = T)
            tot <- cbind(set = "Total", tot)

            df <- pivot_longer(tot, cols = c(4:9))

            p <- df %>% filter(name == "QD") %>%
                    ggplot(aes(x=value, fill=set)) +
                    geom_density(alpha=0.4) +
                    geom_vline(aes(xintercept=2), color = "purple", linetype="dashed") +
                    xlab("QD")

            ggsave(p, filename = "QD.png")

            p <- df %>% filter(name == "FS") %>%
                    ggplot(aes(x=value, fill=set)) +
                        geom_density(alpha=0.4) +
                        xlab("FS") + geom_vline(aes(xintercept=60), color = "purple", linetype="dashed")

            ggsave(p, filename = "FS.png")

            p <- df %>% filter(name == "SOR") %>%
                    ggplot(aes(x=value, fill=set)) +
                        geom_density(alpha=0.4) +
                        xlab("SOR") +
                        geom_vline(aes(xintercept=3), color = "purple", linetype="dashed")

            ggsave(p, filename = "SOR.png")

            p <- df %>% filter(name == "MQ") %>%
                    ggplot(aes(x=value, fill=set)) +
                        geom_density(alpha=0.4) +
                        xlab("MQ") +
                        geom_vline(aes(xintercept=40), color = "purple", linetype="dashed")

            ggsave(p, filename = "MQ.png")

            p <- df %>% filter(name == "MQRankSum") %>%
                    ggplot(aes(x=value, fill=set)) +
                        geom_density(alpha=0.4) +
                        xlab("MQRankSum") +
                        geom_vline(aes(xintercept=-12.5), color = "purple", linetype="dashed")

            ggsave(p, filename = "MQRankSum.png")

            p <- df %>% filter(name == "ReadPosRankSum") %>%
                    ggplot(aes(x=value, fill=set)) +
                        geom_density(alpha=0.4) +
                        xlab("ReadPosRankSum") +
                        geom_vline(aes(xintercept=-8), color = "purple", linetype="dashed")

            ggsave(p, filename = "ReadPosRankSum.png")

            system("mkdir QualPlots")
            system("mv *png QualPlots")
            system("tar -czvf QualPlots.tar.gz QualPlots")

        RSCRIPT

    >>>

    runtime {
        docker: "cristaniguti/reads2map:0.0.4"
        cpu: 1
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "QualPlots"
        mem:"~{memory_size}M"
        time:"01:40:00"
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Generated graphics about empirical markers quality parameters and Hard Filtering. See more information in [VariatsToTable](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants) tool"
    }

    output {
        File Plots = "QualPlots.tar.gz"
    }
}


task QualPlotsForHardFilteringSimulated {
    input {
        File FalsePositives
        File TruePositives
        File Total
        Int? seed
        Int? depth
    }

    Int disk_size = ceil(size(FalsePositives, "GB") + size(TruePositives, "GB") + size(Total, "GB") + 1)
    Int memory_size = ceil(size(Total, "MiB") * 2) + 4000

    command <<<
        R --vanilla --no-save <<RSCRIPT
            system("cp ~{FalsePositives} .")
            system("cp ~{TruePositives} .")
            system("cp ~{Total} .")

            library(ggplot2)
            library(dplyr)
            library(tidyr)

            tot <- read.table("Total.table", header = T)
            tot <- cbind(set = "Total", tot)

            fp <- read.table("FalsePositives.table", header = T)
            fp <- cbind(set = "False positive", fp)

            tp <- read.table("TruePositives.table", header = T)
            tp <- cbind(set = "True positive", tp)

            df <- rbind(tot, fp, tp)
            df <- pivot_longer(df, cols = c(4:9))

            p <- df %>% filter(name == "QD") %>%
                    ggplot(aes(x=value, fill=set)) +
                    geom_density(alpha=0.4) +
                    geom_vline(aes(xintercept=2), color = "purple", linetype="dashed") +
                    xlab("QD")

            ggsave(p, filename = "QD.png")

            p <- df %>% filter(name == "FS") %>%
                    ggplot(aes(x=value, fill=set)) +
                        geom_density(alpha=0.4) +
                        xlab("FS") + geom_vline(aes(xintercept=60), color = "purple", linetype="dashed")

            ggsave(p, filename = "FS.png")

            p <- df %>% filter(name == "SOR") %>%
                    ggplot(aes(x=value, fill=set)) +
                        geom_density(alpha=0.4) +
                        xlab("SOR") +
                        geom_vline(aes(xintercept=3), color = "purple", linetype="dashed")

            ggsave(p, filename = "SOR.png")

            p <- df %>% filter(name == "MQ") %>%
                    ggplot(aes(x=value, fill=set)) +
                        geom_density(alpha=0.4) +
                        xlab("MQ") +
                        geom_vline(aes(xintercept=40), color = "purple", linetype="dashed")

            ggsave(p, filename = "MQ.png")

            p <- df %>% filter(name == "MQRankSum") %>%
                    ggplot(aes(x=value, fill=set)) +
                        geom_density(alpha=0.4) +
                        xlab("MQRankSum") +
                        geom_vline(aes(xintercept=-12.5), color = "purple", linetype="dashed")

            ggsave(p, filename = "MQRankSum.png")

            p <- df %>% filter(name == "ReadPosRankSum") %>%
                    ggplot(aes(x=value, fill=set)) +
                        geom_density(alpha=0.4) +
                        xlab("ReadPosRankSum") +
                        geom_vline(aes(xintercept=-8), color = "purple", linetype="dashed")

            ggsave(p, filename = "ReadPosRankSum.png")

            system("mkdir ~{seed}_~{depth}_QualPlots")
            system("mv *png ~{seed}_~{depth}_QualPlots")
            system("tar -czvf ~{seed}_~{depth}_QualPlots.tar.gz ~{seed}_~{depth}_QualPlots")

        RSCRIPT

    >>>

    runtime {
        docker: "cristaniguti/reads2map:0.0.4"
        cpu: 1
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "QualPlots"
        mem:"~{memory_size}M"
        time:"01:40:00"
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Generated graphics about simulated markers quality parameters and Hard Filtering. See more information in [VariatsToTable](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants) tool"
    }

    output {
        File Plots = "~{seed}_~{depth}_QualPlots.tar.gz"
    }

}

task FilterMulti {
    input {
        File multi_vcf
        String? P1
        String? P2
        Int ploidy
    }

    Int disk_size = ceil(size(multi_vcf, "GiB") * 1.5)
    Int memory_size = 3000

    command <<<
        R --vanilla --no-save <<RSCRIPT

            library(Reads2MapTools)
            filter_multi_vcf("~{multi_vcf}", "~{P1}", "~{P2}",
                             ploidy = ~{ploidy},
                             vcf.out = "multi_vcf_filt.vcf.gz")

        RSCRIPT
    >>>

    runtime {
        docker:"cristaniguti/reads2map:0.0.4"
        cpu: 1
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "FilterMulti"
        mem:"~{memory_size}M"
        time:"01:00:00"
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Filters VCF file markers according to segregation expected in a outcrossing F1 population. Adapts the alleles codification. See [Reads2MapTools](https://github.com/Cristianetaniguti/Reads2MapTools) for more information."
    }

    output {
        File multi_vcf_filt = "multi_vcf_filt.vcf.gz"
    }
}
