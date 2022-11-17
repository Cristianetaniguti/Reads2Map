version 1.0

import "../structs/read_simulation_structs.wdl"
import "../tasks/BWA.wdl" as alg
import "../tasks/chunk_lists.wdl"
import "../tasks/simuscop.wdl"
import "../tasks/pedigree_simulator_utils.wdl"
import "../tasks/pirs.wdl"
import "../tasks/pedigree_simulator.wdl"
import "../tasks/utils.wdl" as utils
import "../tasks/radinitio.wdl"

workflow CreateAlignmentFromSimulation {
    input {
        ReferenceFasta references
        Family family
        Sequencing sequencing
        Int max_cores
        Int chunk_size
    }

  # User can provide specific variants in VCF file
  # If not, use pirs to simulate
  if (!defined(sequencing.emp_vcf)){
    call pirs.GenerateAlternativeGenome {
      input:
        seed       = family.seed,
        ref_genome = references.ref_fasta
    }

    call pedigree_simulator_utils.CreatePedigreeSimulatorInputs {
      input:
        seed       = family.seed,
        snps       = GenerateAlternativeGenome.snps,
        indels     = GenerateAlternativeGenome.indels,
        cmBymb     = family.cmBymb,
        ref        = references.ref_fasta,
        ref_fai    = references.ref_fasta_index,
        cross      = family.cross,
        popsize    = family.popsize,
        ploidy     = family.ploidy,
        doses      = family.doses
    }
  }

  if (defined(sequencing.emp_vcf)){
    call pedigree_simulator_utils.Vcf2PedigreeSimulator {
      input:
        vcf_file = sequencing.emp_vcf,
        ref_map = sequencing.ref_map,
        seed = family.seed,
        popsize = family.popsize,
        vcf_parent1 = sequencing.vcf_parent1,
        vcf_parent2 = sequencing.vcf_parent2
    }
  }

  File mapfile_sele = select_first([Vcf2PedigreeSimulator.mapfile_map, CreatePedigreeSimulatorInputs.mapfile_nomap])
  File founderfile_sele = select_first([Vcf2PedigreeSimulator.founderfile_map, CreatePedigreeSimulatorInputs.founderfile_nomap])
  File chromfile_sele = select_first([Vcf2PedigreeSimulator.chromfile_map, CreatePedigreeSimulatorInputs.chromfile_nomap])
  File parfile_sele = select_first([Vcf2PedigreeSimulator.parfile_map, CreatePedigreeSimulatorInputs.parfile_nomap])
  File ref_alt_alleles_sele = select_first([Vcf2PedigreeSimulator.ref_alt_alleles_map, CreatePedigreeSimulatorInputs.ref_alt_alleles_nomap])
  File simulated_phases_sele = select_first([Vcf2PedigreeSimulator.simulated_phases_map, CreatePedigreeSimulatorInputs.simulated_phases_nomap])

  call pedigree_simulator.RunPedigreeSimulator {
    input:
      mapfile     = mapfile_sele,
      founderfile = founderfile_sele,
      chromfile   = chromfile_sele,
      parfile     = parfile_sele
  }

  call pedigree_simulator_utils.ConvertPedigreeSimulationToVcf {
    input:
      seed            = family.seed,
      depth           = sequencing.depth,
      genotypes_dat   = RunPedigreeSimulator.genotypes_dat,
      map_file        = mapfile_sele,
      chrom_file      = chromfile_sele,
      ref_alt_alleles = ref_alt_alleles_sele,
      popsize         = family.popsize,
      mapsize         = sequencing.mapsize
  }

  call utils.GenerateSampleNames {
    input:
      simulated_vcf = ConvertPedigreeSimulationToVcf.simu_vcf
  }

  if(sequencing.library_type == "WGS" || sequencing.library_type == "exome"){

    call simuscop.SimuscopProfile {
      input:
        library_type = sequencing.library_type,
        emp_bam = sequencing.emp_bam,
        vcf     = ConvertPedigreeSimulationToVcf.simu_vcf,
        references = references
    }

    scatter (sampleName in GenerateSampleNames.names) {

      call simuscop.SimuscopSimulation {
        input:
          library_type  = sequencing.library_type,
          sampleName    = sampleName,
          depth         = sequencing.depth,
          emp_bam       = sequencing.emp_bam,
          vcf           = ConvertPedigreeSimulationToVcf.simu_vcf,
          references    = references,
          chrom         = sequencing.chromosome,
          profile       = SimuscopProfile.profile
      }
    }
  }

    # Two option of RADseq
    # The samples need to be simulated together, otherwise they will be all heterozygous
    if(sequencing.library_type == "sdRAD" || sequencing.library_type == "ddRAD"){
      call radinitio.RADinitioSimulation {
        input:
          depth          = sequencing.depth,
          depth_parents  = sequencing.depth_parents,
          enzyme1        = sequencing.enzyme1,
          enzyme2        = sequencing.enzyme2,
          simu_vcf       = ConvertPedigreeSimulationToVcf.simu_vcf,
          radinitio_vcf  = ConvertPedigreeSimulationToVcf.radinitio_vcf,
          references     = references,
          pcr_cycles     = sequencing.pcr_cycles,
          insert_size    = sequencing.insert_size,
          insert_size_dev = sequencing.insert_size_dev,
          read_length    = sequencing.read_length,
          library_type   = sequencing.library_type,
          chrom          = sequencing.chromosome,
          names          = GenerateSampleNames.names
      }
    }


  Array[File] fastq = select_first([RADinitioSimulation.fastq_rad, SimuscopSimulation.fastq_seq])

  call chunk_lists.SepareChunksFastq {
    input:
      fastqs = fastq,
      chunk_size = chunk_size
  }

  scatter (chunk in SepareChunksFastq.chunks){

    call alg.RunBwaAlignmentSimu {
      input:
        reads      = chunk,
        fastqs     = fastq,
        references = references,
        max_cores  = max_cores,
        rm_dupli   = sequencing.rm_dupli
    }
  }

    # Store for MCHap
    call utils.MergeBams {
        input:
            bam_files = flatten(RunBwaAlignmentSimu.bam)
    }

  output {
      Array[File] bam = flatten(RunBwaAlignmentSimu.bam)
      Array[File] bai = flatten(RunBwaAlignmentSimu.bai)
      File ref_alt_alleles = ref_alt_alleles_sele
      Array[String] names = GenerateSampleNames.names
      File true_vcf = ConvertPedigreeSimulationToVcf.simu_vcf
      File simu_haplo = ConvertPedigreeSimulationToVcf.simu_haplo
      File simulated_phases = simulated_phases_sele
      File merged_bam = MergeBams.merged_bam
  }
}
