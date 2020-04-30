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

struct Dataset {
  File samples_info
  String name
  String parent1
  String parent2
  String? filt_depth
  String cross
}
