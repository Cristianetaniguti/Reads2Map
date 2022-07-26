version 1.0

import "./utils.wdl" as utils

workflow gusmapMaps {
  input {
    File vcf_file
    File new_vcf_file
    String SNPCall_program
    String GenotypeCall_program
    String parent1
    String parent2
    Int max_cores
  }

  Array[String] counts                      = ["vcf", "bam"]
  Array[File] vcfs                          = [vcf_file, new_vcf_file]
  Array[Pair[String, File]] counts_and_vcfs = zip(counts, vcfs)

  scatter (vcf in counts_and_vcfs) {
    call GusmapReport{
        input:
          vcf_file = vcf.right,
          SNPCall_program = SNPCall_program,
          GenotypeCall_program = GenotypeCall_program,
          CountsFrom = vcf.left,
          parent1 = parent1,
          parent2 = parent2,
          max_cores = max_cores
        }
    }

    call utils.CompressGusmap{
      input:
        name = "gusmap_map",
        RDatas = GusmapReport.maps_RData,
        maps_report = GusmapReport.maps_report,
        times = GusmapReport.times
    }

   output{
     File tar_gz_report = CompressGusmap.tar_gz_report
   }
}

task GusmapReport{
  input{
    File vcf_file
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    String parent1
    String parent2
    Int max_cores
  }

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
    # preemptible: 3
    # memory:"8 GB"
    # cpu:4
    job_name: "GusmapReport"
    node:"--nodes=1"
    mem:"--mem=10GB"
    tasks:"--ntasks=1"
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
