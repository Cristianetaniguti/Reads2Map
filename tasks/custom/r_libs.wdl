version 1.0

import "../../structs/struct_reference.wdl"

task SimuscopProfile{
  input {
    String library_type
    File?  emp_bam
    File   vcf
    ReferenceFasta  references
  }

  Int disk_size = ceil(size(vcf, "GiB") * 2 + size(emp_bam, "GiB"))
  Int memory_size = 5000

  command <<<
    R --vanilla --no-save <<RSCRIPT

    library(vcfR)
    library(simuscopR)

    if("~{library_type}" == "exome"){

      system(paste0("bamtobed -i ~{emp_bam} > bed_file"))

      seqToProfile("~{emp_bam}", "bed_file", "~{vcf}",
             "~{references.ref_fasta}", "profile")

    } else {
      seqToProfile("~{emp_bam}", vcf.file =  "~{vcf}",
             reference = "~{references.ref_fasta}",  out.profile = "sample.profile")
    }

    RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/reads2map:0.0.1"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "SimuscopProfile"
    mem:"~{memory_size}M"
    time:"05:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Run [seqToProfile](https://github.com/qasimyu/simuscop) to generate data set profile."
  }

  output{
    File profile = "sample.profile"
  }
}

task SimuscopSimulation{
 input {
    String library_type
    String sampleName
    Int depth
    File? emp_bam
    File vcf
    ReferenceFasta references
    String chrom
    File profile
  }

  Int disk_size = ceil(size(emp_bam, "GiB") + size(vcf, "GiB") + size(references.ref_fasta, "GiB") * depth)
  Int memory_size = 10000

  command <<<
    R --vanilla --no-save <<RSCRIPT
      vcfR.object <- read.vcfR("~{vcf}")

      variants <- vcf2variants(vcfR.object, sample = "~{sampleName}", chrom = "~{chrom}")

      write.table(variants$SNVs, file = "SNVs.txt", sep = "\t", quote = F, col.names = F, row.names = F)
      write.table(variants$indels, file = "indels.txt", sep = "\t", quote = F, col.names = F, row.names = F)
      write.table(variants$insertions, file = "insertions.txt", sep = "\t", quote = F, col.names = F, row.names = F)

      system("cat SNVs.txt indels.txt insertions.txt > variants.txt")

      if("~{library_type}" == "exome"){
        simuReads(ref = "~{references.ref_fasta}",
              profile = "~{profile}",
              variation = "variants.txt",
              target = "bed_file",
              name = "~{sampleName}",
              output = ".",
              layout = "SE",
              threads = 6,
              verbose = 1,
              coverage = ~{depth})
      } else {
        simuReads(ref = "~{references.ref_fasta}",
              profile = "profile",
              variation = "variants.txt",
              name = "~{sampleName}",
              output = ".",
              layout = "SE", # only single-end by now
              threads = 6,
              verbose = 1,
              coverage = ~{depth})
      }

    RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/reads2map:0.0.1"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "SimuscopSimulation"
    mem:"~{memory_size}M"
    time:"10:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Run [simuReads](https://github.com/qasimyu/simuscop) to simulated exome or WGS sequencing reads."
  }

  output {
    File fastq_seq = "~{sampleName}.fq"
  }
}


task GusmapReport {
  input{
    File vcf_file
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    String parent1
    String parent2
    Int max_cores
  }

  Int disk_size = ceil(size(vcf_file, "GiB") * 1.5)
  Int memory_size = ceil(size(vcf_file, "MiB"))

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

  runtime{
    docker:"cristaniguti/reads2map:0.0.1"
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

  output{
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

  Int disk_size = ceil(size(vcf_file, "GiB") * 1.5 + size(simu_onemap_obj, "GiB") + size(ref_alt_alleles, "GiB") + size(simulated_phases, "GiB") + 3)
  Int memory_size = ceil(size(vcf_file, "MiB"))

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
    docker: "cristaniguti/reads2map:0.0.1"
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

task QualPlots {
    input {
        File Total
    }

    Int disk_size = ceil(size(Total, "GB") + 1)
    Int memory_size = ceil(size(Total, "MiB") * 1.25)

    command <<<
        R --vanilla --no-save <<RSCRIPT
            system("cp ~{Total} .")

            library(ggplot2)
            library(dplyr)
            library(tidyr)

            tot <- read.table("Total.table", header = T)
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
        docker: "cristaniguti/reads2map:0.0.1"
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
    Int memory_size = ceil(size(Total, "MiB") * 1.25)

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
        docker: "cristaniguti/reads2map:0.0.1"
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


task VariantFiltration {
    input {
        File vcf_file
        File vcf_tbi
        File reference
        File reference_idx
        File reference_dict
    }

    Int disk_size = ceil(size(vcf_file, "GB") + size(reference, "GB") + 1)
    Int memory_size = 5000

    command <<<
        /usr/gitc/gatk4/./gatk VariantFiltration \
            -V ~{vcf_file} \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O gatk_filters.vcf.gz

        /usr/gitc/gatk4/./gatk SelectVariants \
            -R ~{reference} \
            -V gatk_filters.vcf.gz \
            --exclude-filtered \
            -O gatk_filtered.vcf.gz

    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
        cpu: 1
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "VariantFiltration"
        mem:"~{memory_size}M"
        time:"01:00:00"
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Filters simulated VCF according to GATK Hard Filtering. See more information in [VariatsToTable](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants) tool"
    }

    output {
        File filters_vcf = "gatk_filters.vcf.gz"
        File filtered_vcf = "gatk_filtered.vcf.gz"
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
    Int memory_size = 10000

    command <<<
        R --vanilla --no-save <<RSCRIPT

            library(Reads2MapTools)
            filter_multi_vcf("~{multi_vcf}", "~{P1}", "~{P2}",
                             ploidy = ~{ploidy},
                             vcf.out = "multi_vcf_filt.vcf.gz")

        RSCRIPT
    >>>

    runtime {
        docker:"cristaniguti/reads2map:0.0.1"
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


task JointReports {
  input {
    Array[File] SNPCaller
    Array[File] updog
    Array[File] polyrad
    Array[File] supermassa
    Array[File] gusmap_files
    Array[File] multiallelics_file
    File Freebayes_eval
    File GATK_eval
    Int max_cores
    Int seed
    Int depth
  }

  Int disk_size = ceil(size(SNPCaller, "GiB") * 1.5 + size(updog, "GiB") * 1.5 + size(polyrad, "GiB") * 1.5 + size(supermassa, "GiB") * 1.5 + size(gusmap_files, "GiB") * 1.5 + size(multiallelics_file, "GiB") * 1.5)
  Int memory_size = 5000

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
      gusmap <- str_split("~{sep=";" gusmap_files}", ";", simplify = T)

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
      eval1 <- cbind(SNPCall = "Freebayes", seed = ~{seed}, depth = ~{depth}, df[["CompOverlap"]])
      count1 <- cbind(SNPCall = "Freebayes", seed = ~{seed}, depth = ~{depth}, df[["CountVariants"]])

      df <- gsa.read.gatkreport("~{GATK_eval}")
      eval2 <- cbind(SNPCall = "GATK", seed = ~{seed}, depth = ~{depth}, df[["CompOverlap"]])
      count2 <- cbind(SNPCall = "GATK", seed = ~{seed}, depth = ~{depth}, df[["CountVariants"]])

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
    cpu: max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "JointReports"
    mem:"~{memory_size}M"
    time:"01:40:00"
  }

  meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Merge reports resulted from evaluations with different pipelines."
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
