version 1.0

struct Family {
  Float? cmBymb
  Int seed
  Int popsize
  File? doses
  Int ploidy
  String cross
}

struct Sequencing {
  String library_type
  Int depth
  Int depth_parents
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
  String rm_dupli
  Int mapsize
}
