version 1.0

import "./utilsR.wdl" as utilsR

workflow SimulatedMap{
  input{
    File vcf_simu
    String cross
  }

  call utilsR.vcf2onemap{
    input:
      vcf_file = vcf_simu,
      cross = cross,
      SNPCall_program = "simu"
  }

  output{
    File simu_onemap_obj = vcf2onemap.onemap_obj
  }
}