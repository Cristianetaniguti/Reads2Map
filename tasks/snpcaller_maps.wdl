version 1.0

import "./utilsR.wdl" as utilsR

workflow SNPCallerMaps{
  input {
     File simu_onemap_obj
     File onemap_obj
     File vcfR_obj
     File tot_mks
     File real_phases
     String cross
     String SNPCall_program
     String GenotypeCall_program
     String CountsFrom
     String cMbyMb
    }
    
  call utilsR.GQProbs{
    input:
      vcfR_obj = vcfR_obj,
      onemap_obj = onemap_obj,
      cross = cross
  }
  
  call utilsR.FiltersReport{
    input:
      onemap_obj = GQProbs.gq_onemap_obj,
      SNPCall_program = SNPCall_program,
      GenotypeCall_program = "SNPCaller",
      CountsFrom = "vcf"
  }
            
  call utilsR.MapsReport{
    input:
      onemap_obj = FiltersReport.onemap_obj_filtered,
      tot_mks = tot_mks,
      simu_onemap_obj = simu_onemap_obj,
      SNPCall_program = SNPCall_program,
      GenotypeCall_program = "SNPCaller",
      CountsFrom = CountsFrom,
      cMbyMb = cMbyMb,
      real_phases = real_phases
  }
            
  call utilsR.ErrorsReport{
    input:
      onemap_obj = GQProbs.gq_onemap_obj,
      simu_onemap_obj = simu_onemap_obj,
      SNPCall_program = SNPCall_program,
      GenotypeCall_program = "SNPCaller",
      CountsFrom = CountsFrom
  }
  
  output{
    File RDatas = MapsReport.maps_RData
    File maps_report = MapsReport.maps_report
    File filters_report = FiltersReport.filters_report
    File errors_report = ErrorsReport.errors_report
  }
}