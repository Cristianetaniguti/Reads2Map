version 1.0

struct Samples_info {
  File samples_info
}

struct SplitVCF {
  String chromosome
  String parent1
  String parent2
}

struct References{
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File ref_sa
  File ref_amb
  File ref_bwt
  File ref_ann
  File ref_pac
}