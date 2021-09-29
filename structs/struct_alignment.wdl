version 1.0

struct Alignment {
  File bam
  File bai
  String sample
}


struct Variants {
  File vcf
  File tbi
  String sample
}
