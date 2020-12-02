version 1.0

struct Family {
  String cmBymb
  Int seed
  Int popsize
  File doses
  Int ploidy
  String cross
}

struct Sequencing{
  String library_type
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
}

struct OptionalFilters{
  String? Filter1
  String? Filter2
  String? Filter3
  String? Filter4
  String? Filter5
}
