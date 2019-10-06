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
  String name
  Int popsize
  String enzyme
  Int depth
  File doses
  Int ploidy
  String cross
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
