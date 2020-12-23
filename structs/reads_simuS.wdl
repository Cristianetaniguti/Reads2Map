version 1.0

struct Family {
  Float? cmBymb
  Int seed
  Int popsize
<<<<<<< HEAD
  File? doses
=======
  String enzyme1
  String enzyme2
  Int depth
  File doses
>>>>>>> master
  Int ploidy
  String cross
}

<<<<<<< HEAD
struct Sequencing{
  String library_type
=======
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
>>>>>>> master
  Int depth
  Int? insert_size
  Int? insert_size_dev
  Int? pcr_cycles
  Int? read_length
  String enzyme1
  String? enzyme2
  File emp_vcf
  File? emp_bam
  File? ref_map
  String chromosome
  String multiallelics
  String vcf_parent1
  String vcf_parent2
}

struct OptionalFilters{
  String? Filter1
  String? Filter2
  String? Filter3
  String? Filter4
  String? Filter5
}
