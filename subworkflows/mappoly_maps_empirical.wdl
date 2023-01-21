version 1.0

import "../tasks/utilsR.wdl" as utilsR
import "../tasks/mappoly.wdl" as mappolyTasks

workflow MappolyMapsEmp {

  input {
    File vcf_file
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    String cross
    String parent1
    String parent2
    Int max_cores
    Int ploidy
  }

  call utilsR.ReGenotyping {
      input:
          vcf_file = vcf_file,
          GenotypeCall_program = GenotypeCall_program,
          cross = cross,
          parent1 = parent1,
          parent2 = parent2,
          max_cores = max_cores,
          ploidy = ploidy
  }

  call mappolyTasks.MappolyReport {
    input:
      vcf_file = ReGenotyping.regeno_vcf,
      parent1 = parent1,
      parent2 = parent2,
      GenotypeCall_program = GenotypeCall_program,
      SNPCall_program = SNPCall_program,
      CountsFrom = CountsFrom,
      max_cores = max_cores,
      ploidy = ploidy
  }

   output {
      File tar_gz_report = MappolyReport.results
      File regeno_vcf = ReGenotyping.regeno_vcf
   }
}
