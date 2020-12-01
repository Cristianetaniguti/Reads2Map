version 1.0

struct Family {
  String cmBymb
  Int seed
  Int popsize
  String enzyme1
  String enzyme2
  Int depth
  File doses
  Int ploidy
  String cross
  String multiallelics
  String parent1
  String parent2
}

struct Sequencing{
  String type
  Int depth
  Int? insert_size
  Int? isert_size_dev
  Int? pcr_cycles
  Int? read_length
  String enzyme1
  String? enzyme2
  File emp_vcf
  File? emp_bam
  File? ref_map
  String chromosome
}


struct Profiles{
  File base_calling
  File indel_error
  File gc_bias
}

struct FamilyTemplate {
  String cmBymb
  Int? seed
  Int popsize
  String enzyme1
  String enzyme2
  Int depth
  File doses
  Int ploidy
  String cross
  String chromosome
  Int global_seed
  String multiallelics
}

struct SplitVCF{
  String chromosome
  String parent1
  String parent2
}

struct OptionalFilters{
  String? Filter1
  String? Filter2
  String? Filter3
  String? Filter4
  String? Filter5
}
