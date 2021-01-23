version 1.0

import "simulated_map.wdl" as simulated_map
import "utilsR.wdl" as utilsR
import "default_maps.wdl" as default
import "snpcaller_maps.wdl" as snpcaller
import "genotyping-simulated.wdl" as genotyping
import "gusmap_maps.wdl" as gusmap


workflow TestGenotypingR {

    input {
        File analysis_bam
        File analysis_vcf
        File analysis_multi_vcf
        File true_vcf
        File ref_alt_alleles
        File simulated_phases
        String method
        String parent1
        String parent2
        String cross
        Int seed
        Int depth
        Int max_cores
    }

    call simulated_map.SimulatedMap {
        input:
        vcf_simu = true_vcf,
        cross = cross
    }

    call utilsR.vcf2onemap {
      input:
        vcf_file = analysis_vcf,
        cross = cross,
        SNPCall_program = method,
        parent1 = parent1,
        parent2 = parent2
    }

    call utilsR.MultiVcf2onemap{
       input:
          multi = analysis_multi_vcf,
          cross = cross,
          SNPCall_program = method,
          parent1 = parent1,
          parent2 = parent2,
          seed    = seed,
          depth   = depth
    }

    call default.DefaultMaps {
      input:
        onemap_obj = vcf2onemap.onemap_obj,
        simu_onemap_obj = SimulatedMap.simu_onemap_obj,
        ref_alt_alleles = ref_alt_alleles,
        simulated_phases = simulated_phases,
        SNPCall_program = method,
        CountsFrom = "vcf",
        multi_obj = MultiVcf2onemap.onemap_obj
    }

    call snpcaller.SNPCallerMaps {
      input:
        simu_onemap_obj = SimulatedMap.simu_onemap_obj,
        onemap_obj = vcf2onemap.onemap_obj,
        vcf_file = analysis_vcf,
        ref_alt_alleles = ref_alt_alleles,
        simulated_phases = simulated_phases,
        cross = cross,
        SNPCall_program = method,
        GenotypeCall_program = "SNPCaller",
        CountsFrom = "vcf",
        multi_obj = MultiVcf2onemap.onemap_obj
    }

    Map[String, File] vcfs = {"vcf": analysis_vcf, "bam": analysis_bam}

    scatter (origin in ["vcf", "bam"]){
        call genotyping.SnpBasedGenotypingSimulatedMaps as UpdogMaps {
          input:
            simu_onemap_obj = SimulatedMap.simu_onemap_obj,
            onemap_obj = vcf2onemap.onemap_obj,
            vcf_file = vcfs[origin],
            genotyping_program = "updog",
            ref_alt_alleles = ref_alt_alleles,
            simulated_phases = simulated_phases,
            SNPCall_program = method,
            CountsFrom = origin,
            cross = cross,
            multi_obj = MultiVcf2onemap.onemap_obj,
            max_cores = max_cores
        }

        call genotyping.SnpBasedGenotypingSimulatedMaps as SupermassaMaps {
          input:
            simu_onemap_obj = SimulatedMap.simu_onemap_obj,
            onemap_obj = vcf2onemap.onemap_obj,
            vcf_file = vcfs[origin],
            genotyping_program = "supermassa",
            ref_alt_alleles = ref_alt_alleles,
            simulated_phases = simulated_phases,
            SNPCall_program = method,
            CountsFrom = origin,
            cross = cross,
            multi_obj = MultiVcf2onemap.onemap_obj,
            max_cores = max_cores
        }

        call genotyping.SnpBasedGenotypingSimulatedMaps as PolyradMaps {
          input:
            simu_onemap_obj = SimulatedMap.simu_onemap_obj,
            onemap_obj = vcf2onemap.onemap_obj,
            vcf_file = vcfs[origin],
            genotyping_program = "polyrad",
            ref_alt_alleles = ref_alt_alleles,
            simulated_phases = simulated_phases,
            SNPCall_program = method,
            CountsFrom = origin,
            cross = cross,
            multi_obj = MultiVcf2onemap.onemap_obj,
            max_cores = max_cores
        }
      }

      call gusmap.GusmapMaps {
        input:
          simu_onemap_obj = SimulatedMap.simu_onemap_obj,
          vcf_file = analysis_vcf,
          new_vcf_file = analysis_bam,
          SNPCall_program = method,
          GenotypeCall_program = "gusmap",
          ref_alt_alleles = ref_alt_alleles,
          simulated_phases = simulated_phases
      }

      output {
          Array[File] gusmap_out = GusmapMaps.RDatas
      }
}
