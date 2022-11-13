version 1.0

import "../tasks/utilsR.wdl" as utilsR
import "../tasks/utils.wdl" as utils

workflow SNPCallerMaps {
  input {
     File simu_onemap_obj
     File simu_vcfR
     File vcf_file
     File ref_alt_alleles
     File simulated_phases
     String cross
     String SNPCall_program
     String GenotypeCall_program
     String CountsFrom
     Int seed
     Int depth
     Int max_cores
     String multiallelics
     File? multiallelics_file
     File? multiallelics_mchap
     String mchap
    }

  if (multiallelics == "TRUE") {
    call utils.JointMarkers {
      input:
        biallelic_vcf = vcf_file,
        multiallelic_vcf = multiallelics_file
    }
  }

  File updated_vcf = select_first([JointMarkers.merged_vcf, vcf_file])

  call utilsR.SetProbsDefault {
    input:
      vcf_file = updated_vcf,
      cross = cross,
      parent1 = "P1",
      parent2 = "P2",
      SNPCall_program = SNPCall_program,
      multiallelics = multiallelics,
      multiallelics_mchap = multiallelics_mchap,
      mchap = mchap
  }

  Array[String] methods                         = [GenotypeCall_program, GenotypeCall_program + "0.05", GenotypeCall_program + "default"]
  Array[File] objects                           = [SetProbsDefault.probs_onemap_obj, SetProbsDefault.globalerror_onemap_obj, SetProbsDefault.default_onemap_obj]
  Array[Pair[String, File]] methods_and_objects = zip(methods, objects)

  scatter (item in methods_and_objects) {
      call utilsR.FiltersReport {
        input:
          onemap_obj = item.right,
          SNPCall_program = SNPCall_program,
          GenotypeCall_program = item.left,
          CountsFrom = "vcf",
          seed = seed,
          depth = depth
      }

      call utilsR.MapsReport {
        input:
          onemap_obj = FiltersReport.onemap_obj_filtered,
          ref_alt_alleles = ref_alt_alleles,
          simu_onemap_obj = simu_onemap_obj,
          SNPCall_program = SNPCall_program,
          GenotypeCall_program = item.left,
          CountsFrom = CountsFrom,
          simulated_phases = simulated_phases,
          seed = seed,
          depth = depth,
          max_cores = max_cores
      }

      call utilsR.ErrorsReport {
        input:
          onemap_obj = item.right,
          SNPCall_program = SNPCall_program,
          GenotypeCall_program = item.left,
          CountsFrom = CountsFrom,
          simu_vcfR = simu_vcfR,
          vcfR_obj = SetProbsDefault.vcfR_obj,
          seed = seed,
          depth = depth,
          max_cores = max_cores
      }
  }

  call utils.Compress {
      input:
        RDatas = MapsReport.maps_RData,
        maps_report = MapsReport.maps_report,
        times = MapsReport.times,
        filters_report = FiltersReport.filters_report,
        errors_report = ErrorsReport.errors_report,
        name = "snpcaller_maps"
   }

   output {
      File tar_gz_report = Compress.tar_gz_report
   }
}
