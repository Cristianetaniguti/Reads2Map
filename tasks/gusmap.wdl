version 1.0

task GusmapReport {
  input {
    File vcf_file
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    String parent1
    String parent2
    Int max_cores
  }

  Int disk_size = ceil(size(vcf_file, "GiB") * 2)
  Int memory_size = ceil(size(vcf_file, "MiB") + 2000)

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)
      library(GUSMap)
      library(Reads2MapTools)

      if(tail(strsplit("~{vcf_file}", "[.]")[[1]],1) =="gz") {
          vcf.temp <- paste0("~{SNPCall_program}",".vcf")
          system(paste0("zcat ", "~{vcf_file}", " > ", vcf.temp))
          vcf_file <- vcf.temp
       } else {
          vcf_file <- "~{vcf_file}"
       }

      times_temp <- system.time(info <- create_gusmap_report_emp(vcf_file,"~{SNPCall_program}", "~{CountsFrom}",
                           "~{GenotypeCall_program}", "~{parent1}", "~{parent2}"))

      times <- data.frame(SNPCall = "~{SNPCall_program}",
                    CountsFrom = "~{CountsFrom}",
                    GenoCall =  "~{GenotypeCall_program}",
                    time = times_temp[3])

      vroom::vroom_write(info[[2]], "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_map_report.tsv.gz", num_threads = ~{max_cores})
      vroom::vroom_write(times, "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_times_report.tsv.gz", num_threads = ~{max_cores})
      map_out <- info[[1]]
      save(map_out, file= "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.RData")

    RSCRIPT

  >>>

  runtime {
    docker:"cristaniguti/reads2map:0.0.4"
    cpu: max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "GusmapReport"
    mem:"~{memory_size}G"
    time:"24:00:00"
  }

  meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Estimate genetic distances by GUSMap HMM multi-point approach in a set o markers ordered by genomic position. See [Reads2MapTools](https://github.com/Cristianetaniguti/Reads2MapTools) for more information."
  }

  output {
    File maps_report = "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_map_report.tsv.gz"
    File maps_RData = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.RData"
    File times = "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_times_report.tsv.gz"
  }
}

task GusmapReportForSimulated {
  input {
    File vcf_file
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    File simu_onemap_obj
    File ref_alt_alleles
    File simulated_phases
    Int seed
    Int depth
    Int max_cores
  }

  Int disk_size = ceil(size(vcf_file, "GiB") * 2 + size(simu_onemap_obj, "GiB") + size(ref_alt_alleles, "GiB") + size(simulated_phases, "GiB") + 3)
  Int memory_size = ceil(size(vcf_file, "MiB") + 2000)

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)
      library(GUSMap)
      library(Reads2MapTools)

      simu_onemap_obj <- load("~{simu_onemap_obj}")
      simu_onemap_obj <- get(simu_onemap_obj)

      if(tail(strsplit("~{vcf_file}", "[.]")[[1]],1) =="gz") {
          vcf.temp <- paste0("~{SNPCall_program}",".vcf")
          system(paste0("zcat ", "~{vcf_file}", " > ", vcf.temp))
          vcf_file <- vcf.temp
       } else {
          vcf_file <- "~{vcf_file}"
       }

      ref_alt_alleles <- read.table("~{ref_alt_alleles}")
      simulated_phases <- read.table("~{simulated_phases}")

      times_fake <- system.time(info_fake <- create_gusmap_report_simu(vcf_file, gab= simu_onemap_obj,"~{SNPCall_program}",
                                                     "~{GenotypeCall_program}", fake = "with-false", "~{CountsFrom}", ref_alt_alleles,simulated_phases,
                                                     ~{seed}, ~{depth}))

      times <- data.frame(seed = ~{seed}, depth = ~{depth}, SNPCall = "~{SNPCall_program}",
                          CountsFrom = "~{CountsFrom}", GenoCall =  "~{GenotypeCall_program}", fake = "with-false",
                          time = times_fake[3])

      # If there is no false positive, map will not run again
      if(all(info_fake[[2]][,"real.mks"] == "true marker")){
        cat("skip :) \n")
        times_temp <- times_fake
        info_correct <- update_fake_info(info_fake, simu_onemap_obj, ref_alt_alleles, simulated_phases)

      } else {
        times_temp <- system.time(info_correct <- create_gusmap_report_simu(vcf_file, gab= simu_onemap_obj, "~{SNPCall_program}",
                                                      "~{GenotypeCall_program}", fake = "without-false", "vcf", ref_alt_alleles,simulated_phases,
                                                     ~{seed}, ~{depth}))

      }

      # Joint maps data.frames
      map_joint <- rbind(info_fake[[2]], info_correct[[2]])
      vroom::vroom_write(map_joint, "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_~{seed}_~{depth}_map_report.tsv.gz", num_threads = ~{max_cores})

      # Joint RDatas
      RDatas_joint <- list()
      RDatas_joint[[1]] <- info_fake[[1]]
      RDatas_joint[[2]] <- info_correct[[1]]
      names(RDatas_joint) <- c("map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_TRUE",
                               "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_FALSE")
      save(RDatas_joint, file= "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_~{seed}_~{depth}.RData")

      # Joint times data.frames
      times_temp <- data.frame(seed = ~{seed}, depth = ~{depth}, SNPCall = "~{SNPCall_program}",
                               CountsFrom = "~{CountsFrom}", GenoCall =  "~{GenotypeCall_program}", fake = "without-false",
                               time = times_temp[3])

      times <- rbind(times, times_temp)
      vroom::vroom_write(times, "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_~{seed}_~{depth}_times_report.tsv.gz", num_threads = ~{max_cores})

    RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/reads2map:0.0.4"
    cpu: max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "GusmapReport"
    mem:"~{memory_size}M"
    time:"24:00:00"
  }

  meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Estimate genetic distances by GUSMap HMM multi-point approach in a set o simulated markers ordered by genomic position with and without false-positives. See [Reads2MapTools](https://github.com/Cristianetaniguti/Reads2MapTools) for more information."
  }

  output {
    File maps_report = "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_~{seed}_~{depth}_map_report.tsv.gz"
    File maps_RData = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_~{seed}_~{depth}.RData"
    File times = "~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}_~{seed}_~{depth}_times_report.tsv.gz"
  }
}

task CompressGusmap {
    input{
      String name
      File RDatas
      File maps_report
      File times
    }

    Int disk_size = ceil(size(RDatas, "GiB") + size(maps_report, "GiB") + size(times, "GiB"))
    Int memory_size = 1000

    command <<<

      mkdir ~{name}
      cp ~{RDatas} ~{maps_report} \
                ~{times} ~{name}

      tar -czvf ~{name}.tar.gz ~{name}

    >>>

  runtime {
    docker:"ubuntu:20.04"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "CompressGusmap"
    mem:"~{memory_size}M"
    time:"01:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Move GUSMap resulted reports to a single directory and compress it."
  }

  output {
    File tar_gz_report = "~{name}.tar.gz"
  }
}

task CompressGusmapSimu {
    input{
      String name
      Array[File] RDatas
      Array[File] maps_report
      Array[File] times
    }

    Int disk_size = ceil(size(RDatas, "GiB") + size(maps_report, "GiB") + size(times, "GiB"))
    Int memory_size = 1000

    command <<<

      mkdir ~{name}
      cp ~{sep=" " RDatas} ~{sep=" " maps_report} \
                ~{sep=" " times} ~{name}

      tar -czvf ~{name}.tar.gz ~{name}

    >>>

  runtime {
    docker:"ubuntu:20.04"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "CompressGusmap"
    mem:"~{memory_size}M"
    time:"01:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Move GUSMap resulted reports to a single directory and compress it."
  }

  output {
    File tar_gz_report = "~{name}.tar.gz"
  }
}