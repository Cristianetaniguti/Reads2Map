version 1.0

import "./utils.wdl" as utils

workflow GusmapMaps{
  input{
    File vcf_file
    File new_vcf_file
    String SNPCall_program
    String GenotypeCall_program
    String parent1
    String parent2
    String chromosome
  }

  Array[String] counts                      = ["vcf", "bam"]
  Array[File] vcfs                          = [vcf_file, new_vcf_file]
  Array[Pair[String, File]] counts_and_vcfs = zip(counts, vcfs)
  
  scatter(vcf in counts_and_vcfs){
  
    call utils.SelectChrVCF{
      input:
        vcf_file = vcf.right,
        chromosome = chromosome
    }
    
    call GusmapReport{
        input:
          vcf_file = SelectChrVCF.chr_filt,
          SNPCall_program = SNPCall_program,
          GenotypeCall_program = GenotypeCall_program,
          CountsFrom = vcf.left,
          parent1 = parent1,
          parent2 = parent2
        }
    }
     
   output{
      Array[File] RDatas = GusmapReport.maps_RData
      Array[File] maps_report = GusmapReport.maps_report
      Array[File] times = GusmapReport.times
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
  }
  
  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(onemap)
      library(GUSMap)
      source("/opt/scripts/functions_empirical.R")
      
      if(tail(strsplit("~{vcf_file}", "[.]")[[1]],1) =="gz") {
          vcf.temp <- paste0("~{SNPCall_program}",".vcf")
          system(paste0("zcat ", "~{vcf_file}", " > ", vcf.temp))
          vcf_file <- vcf.temp
       } else {
          vcf_file <- "~{vcf_file}"
       }
      
      create_gusmap_report(vcf_file,"~{SNPCall_program}", "~{CountsFrom}", 
                           "~{GenotypeCall_program}", "~{parent1}", "~{parent2}")

    RSCRIPT
    
  >>>
  
  runtime{
    docker: "taniguti/onemap"
    time:"96:00:00"
    mem:"--nodes=1"
    cpu:1
  }
  
  output{
    File maps_report = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt" 
    File maps_RData = "map_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.RData"
    File times = "times_~{SNPCall_program}_~{CountsFrom}_~{GenotypeCall_program}.txt"
  }
}
