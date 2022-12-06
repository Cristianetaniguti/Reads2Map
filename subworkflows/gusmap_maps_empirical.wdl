version 1.0

import "../tasks/gusmap.wdl" 

workflow gusmapMapsEmp {
  input {
    File vcf_file
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    String parent1
    String parent2
    Int max_cores
  }

    call gusmap.GusmapReport {
        input:
          vcf_file = vcf_file,
          SNPCall_program = SNPCall_program,
          GenotypeCall_program = GenotypeCall_program,
          CountsFrom = CountsFrom,
          parent1 = parent1,
          parent2 = parent2,
          max_cores = max_cores
    }
  
    call gusmap.CompressGusmap {
      input:
        name = "gusmap_map",
        RDatas = GusmapReport.maps_RData,
        maps_report = GusmapReport.maps_report,
        times = GusmapReport.times
    }

   output {
     File tar_gz_report = CompressGusmap.tar_gz_report
   }
}
