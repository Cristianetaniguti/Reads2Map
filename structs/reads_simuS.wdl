version 1.0

struct ReferenceFasta {
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File ref_sa
  File ref_amb
  File ref_bwt
  File ref_ann
  File ref_pac
}

struct Family {
  Float cmBymb
  Int seed
  Int popsize
  String enzyme
  Int depth
  File doses
  Int ploidy
  String cross
}

struct Profiles{
  File base_calling
  File indel_error
  File gc_bias
}

struct FamilyTemplate {
  Float cmBymb
  Int? seed
  Int popsize
  String enzyme
  Int depth
  File doses
  Int ploidy
  String cross
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
