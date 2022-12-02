version 1.0

import "../tasks/gusmap.wdl"

workflow gusmapMaps {
  input {
    File simu_onemap_obj
    File vcf_file
    File new_vcf_file
    String SNPCall_program
    String GenotypeCall_program
    File ref_alt_alleles
    File simulated_phases
    Int seed
    Int depth
    Int max_cores
  }

  Array[String] counts = ["vcf", "bam"]
  Array[File] vcfs = [vcf_file, new_vcf_file]
  Array[Pair[String, File]] counts_and_vcfs = zip(counts, vcfs)

  scatter (vcf in counts_and_vcfs) {
    call gusmap.GusmapReportForSimulated as GusmapReport {
        input:
          vcf_file = vcf.right,
          simu_onemap_obj = simu_onemap_obj,
          SNPCall_program = SNPCall_program,
          GenotypeCall_program = GenotypeCall_program,
          CountsFrom = vcf.left,
          ref_alt_alleles = ref_alt_alleles,
          simulated_phases = simulated_phases,
          seed = seed,
          depth = depth,
          max_cores = max_cores
     }
  }

  call gusmap.CompressGusmapSimu {
     input:
       name = "gusmap_map",
       RDatas = GusmapReport.maps_RData,
       maps_report = GusmapReport.maps_report,
       times = GusmapReport.times
  }

   output {
     File tar_gz_report = CompressGusmapSimu.tar_gz_report
   }
}
