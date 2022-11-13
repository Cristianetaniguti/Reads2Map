version 1.0

import "tasks/custom/pedigree_simulator_utils.wdl"
import "tasks/custom/reports.wdl"

import "subworkflows/SimulatedSingleFamily.wdl" as sub

workflow SimulatedReads {

  input {
    ReferenceFasta references
    Family family
    Sequencing sequencing
    Int number_of_families
    Int global_seed
    Int max_cores
    String? filters
  }

  # ProduceFamiliesSeeds just generates random seeds. It returns an
  # array of integers
  call pedigree_simulator_utils.ProduceFamiliesSeeds {
    input:
      global_seed= global_seed,
      number_of_families=number_of_families
  }

  # Here we generate Family objects on the fly, based on the values
  # from the family and the random seed of the previous task.
  scatter (seed in ProduceFamiliesSeeds.seeds) {
    Family fam =  {
      "cmBymb": family.cmBymb,
      "popsize": family.popsize,
      "enzyme1": sequencing.enzyme1,
      "enzyme2": sequencing.enzyme2,
      "seed": seed,
      "depth": sequencing.depth,
      "doses": family.doses,
      "ploidy": family.ploidy,
      "cross": family.cross,
      "multiallelics": sequencing.multiallelics,
    }

    # Calling reads_simu for each seed
    call sub.SimulatedSingleFamily {
      input:
        references=references,
        family=fam,
        sequencing = sequencing,
        max_cores = max_cores,
        filters = filters,
        ploidy =  family.ploidy
    }
  }

  call reports.JointTables {
    input:
      data1_depths_geno_prob   = SimulatedSingleFamily.data1_depths_geno_prob,
      data2_maps               = SimulatedSingleFamily.data2_maps,
      data3_filters            = SimulatedSingleFamily.data3_filters,
      data5_SNPCall_efficiency = SimulatedSingleFamily.data5_SNPCall_efficiency,
      data4_times              = SimulatedSingleFamily.data4_times,
      data6_RDatas             = SimulatedSingleFamily.data6_RDatas,
      data7_gusmap             = SimulatedSingleFamily.data7_gusmap,
      data8_names              = SimulatedSingleFamily.data8_names,
      data9_simu_haplo         = SimulatedSingleFamily.simu_haplo,
      data10_counts            = SimulatedSingleFamily.data10_counts,
      depth                    = sequencing.depth,
      plots                    = SimulatedSingleFamily.Plots,
      positions                = SimulatedSingleFamily.positions
  }

  # Here you can reference outputs from the sub workflow. Remember that
  # it will be an array of the same type of the original.
  output {
    File results = JointTables.results
  }
}
