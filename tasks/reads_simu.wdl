version 1.0

import "../structs/reads_simuS.wdl"
import "./create_alignment_from_read_simulations.wdl" as simulation
import "./gatk_genotyping.wdl" as gatk
import "./freebayes_genotyping.wdl" as freebayes
import "./utils.wdl" as utils
import "./utilsR.wdl" as utilsR
import "./simulated_map.wdl" as simulated_map
import "./default_maps.wdl" as default
import "./snpcaller_maps.wdl" as snpcaller
import "./updog_maps.wdl" as updog
import "./polyrad_maps.wdl" as polyrad
import "./supermassa_maps.wdl" as supermassa
import "./updogbam_maps.wdl" as updogbam
import "./polyradbam_maps.wdl" as polyradbam
import "./supermassabam_maps.wdl" as supermassabam


workflow reads_simu {

  input {
    ReferenceFasta references
    Family family
    Profiles profiles
  }

  call simulation.CreateAlignmentFromSimulation {
    input:
      references=references,
      family=family,
      profiles=profiles
  }

  call gatk.GatkGenotyping {
    input:
      alignments=CreateAlignmentFromSimulation.alignments,
      references=references,
      program="gatk"
  }

  call freebayes.FreebayesGenotyping {
    input:
      alignments=CreateAlignmentFromSimulation.alignments,
      bam=CreateAlignmentFromSimulation.bam,
      bai=CreateAlignmentFromSimulation.bai,
      references=references,
      program="freebayes"
  }

  call utils.CalculateVcfMetrics {
    input:
      freebayesVCF  = FreebayesGenotyping.vcf,
      gatkVCF       = GatkGenotyping.vcf,
      tot_mks       = CreateAlignmentFromSimulation.total_markers,
      maternal_trim = CreateAlignmentFromSimulation.maternal_trim,
      seed          = family.seed,
      depth         = family.depth
  }

  call utils.BamCounts4Onemap {
    input:
      sampleName       = CreateAlignmentFromSimulation.names,
      freebayes_counts = FreebayesGenotyping.counts,
      gatk_counts      = GatkGenotyping.counts
  }
  
  call simulated_map.SimulatedMap{
    input:
      vcf_simu = CreateAlignmentFromSimulation.true_vcf,
      cross = family.cross
  }

  Array[String] methods                     = ["gatk", "freebayes"]
  Array[File] vcfs                          = [GatkGenotyping.vcf, FreebayesGenotyping.vcf]
  Array[Pair[String, File]] program_and_vcf = zip(methods, vcfs)

  scatter (vcf in program_and_vcf){
    call utilsR.vcf2onemap{
      input:
        vcf_file = vcf.right,
        cross = family.cross,
        SNPCall_program = vcf.left,
        GenotypeCall_program = "default"
    }
    
    call default.DefaultMaps{
      input:
        simu_onemap_obj = SimulatedMap.simu_onemap_obj,
        onemap_obj = vcf2onemap.onemap_obj,
        tot_mks = CreateAlignmentFromSimulation.total_markers,
        real_phases = CreateAlignmentFromSimulation.real_phases,
        SNPCall_program = vcf.left,
        GenotypeCall_program = "default",
        CountsFrom = "vcf",
        cMbyMb = family.cmBymb
    }
    
    call snpcaller.SNPCallerMaps{
      input:
        simu_onemap_obj = SimulatedMap.simu_onemap_obj,
        vcfR_obj = vcf2onemap.vcfR_obj,
        onemap_obj = vcf2onemap.onemap_obj,
        tot_mks = CreateAlignmentFromSimulation.total_markers,
        real_phases = CreateAlignmentFromSimulation.real_phases,
        cross = family.cross,
        SNPCall_program = vcf.left,
        GenotypeCall_program = "SNPCaller",
        CountsFrom = "vcf",
        cMbyMb = family.cmBymb
    }

      call updog.UpdogMaps{
        input:
          simu_onemap_obj = SimulatedMap.simu_onemap_obj,
          vcfR_obj = vcf2onemap.vcfR_obj,
          onemap_obj = vcf2onemap.onemap_obj,
          tot_mks = CreateAlignmentFromSimulation.total_markers,
          real_phases = CreateAlignmentFromSimulation.real_phases,
          SNPCall_program = vcf.left,
          GenotypeCall_program = "updog",
          CountsFrom = "vcf",
          cMbyMb = family.cmBymb
      }
      
      call supermassa.SupermassaMaps{
          input:
          simu_onemap_obj = SimulatedMap.simu_onemap_obj,
          vcfR_obj = vcf2onemap.vcfR_obj,
          onemap_obj = vcf2onemap.onemap_obj,
          tot_mks = CreateAlignmentFromSimulation.total_markers,
          real_phases = CreateAlignmentFromSimulation.real_phases,
          SNPCall_program = vcf.left,
          GenotypeCall_program = "supermassa",
          CountsFrom = "vcf",
          cMbyMb = family.cmBymb
      }
      
      call polyrad.PolyradMaps{
         input:
          simu_onemap_obj = SimulatedMap.simu_onemap_obj,
          vcf_file = vcf.right,
          onemap_obj = vcf2onemap.onemap_obj,
          tot_mks = CreateAlignmentFromSimulation.total_markers,
          real_phases = CreateAlignmentFromSimulation.real_phases,
          SNPCall_program = vcf.left,
          GenotypeCall_program = "polyrad",
          CountsFrom = "vcf",
          cMbyMb = family.cmBymb,
          cross = family.cross
      }
      
      call updogbam.UpdogBamMaps{
        input:
          simu_onemap_obj = SimulatedMap.simu_onemap_obj,
          vcfR_obj = vcf2onemap.vcfR_obj,
          onemap_obj = vcf2onemap.onemap_obj,
          tot_mks = CreateAlignmentFromSimulation.total_markers,
          real_phases = CreateAlignmentFromSimulation.real_phases,
          SNPCall_program = vcf.left,
          GenotypeCall_program = "updog_bam",
          CountsFrom = "bam",
          cMbyMb = family.cmBymb,
          freebayes_ref_bam = BamCounts4Onemap.freebayes_ref_bam,
          freebayes_alt_bam = BamCounts4Onemap.freebayes_alt_bam,
          gatk_ref_bam = BamCounts4Onemap.gatk_ref_bam,
          gatk_alt_bam = BamCounts4Onemap.gatk_alt_bam
      }
      
      call supermassabam.SupermassaBamMaps{
          input:
          simu_onemap_obj = SimulatedMap.simu_onemap_obj,
          vcfR_obj = vcf2onemap.vcfR_obj,
          onemap_obj = vcf2onemap.onemap_obj,
          tot_mks = CreateAlignmentFromSimulation.total_markers,
          real_phases = CreateAlignmentFromSimulation.real_phases,
          SNPCall_program = vcf.left,
          GenotypeCall_program = "supermassa_bam",
          CountsFrom = "bam",
          cMbyMb = family.cmBymb,
          freebayes_ref_bam = BamCounts4Onemap.freebayes_ref_bam,
          freebayes_alt_bam = BamCounts4Onemap.freebayes_alt_bam,
          gatk_ref_bam = BamCounts4Onemap.gatk_ref_bam,
          gatk_alt_bam = BamCounts4Onemap.gatk_alt_bam
      }
      
      call polyradbam.PolyradBamMaps{
         input:
          simu_onemap_obj = SimulatedMap.simu_onemap_obj,
          vcf_file = vcf.right,
          onemap_obj = vcf2onemap.onemap_obj,
          tot_mks = CreateAlignmentFromSimulation.total_markers,
          real_phases = CreateAlignmentFromSimulation.real_phases,
          SNPCall_program = vcf.left,
          GenotypeCall_program = "polyrad_bam",
          CountsFrom = "bam",
          cMbyMb = family.cmBymb,
          cross = family.cross,
          freebayes_ref_bam = BamCounts4Onemap.freebayes_ref_bam,
          freebayes_alt_bam = BamCounts4Onemap.freebayes_alt_bam,
          gatk_ref_bam = BamCounts4Onemap.gatk_ref_bam,
          gatk_alt_bam = BamCounts4Onemap.gatk_alt_bam
      }
  }
}
