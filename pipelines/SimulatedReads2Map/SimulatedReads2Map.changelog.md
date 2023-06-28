# 1.0.2

* more flexibility to choose the probability to be used in the HMM: 

* new parameters:
- global_errors: array with global errors to be tested
- genoprob_error: boolean defining if software probabilities should be tested
- genoprob_global_errors: array with global errors to be tested together with the software probabilities following: 1 - (1 - global error) x (1 - software error probability)

# 1.0.1

* runtimes adapted to run with Caper

# 1.0.0

Initial release

This workflow perform simulations of one or more (defined by `number_of_families`) bi-parental outcrossing population haplotypes using PedigreeSim software based on a provided linkage map and SNP markers. It uses RADinitio software, the simulated haplotypes and a reference genome to also simulate genotyping-by-sequencing read sequences. After, it performs the SNP and genotype calling and builds 68 linkage maps from the combinations:

* SNP calling: GATK and Freebayes
* Dosage/genotype calling: updog, polyRAD and SuperMASSA
* Linkage map build software: OneMap 3.0 and GUSMap
* Using genotype probabilities from GATK, Freebayes, updog, polyRAD and SuperMASSA, and a global error rate of 5% and 0.001% in the OneMap HMM.

It also has the options to:

* Include or not multiallelic (MNP) markers
* Apply filters using VCFtools

This workflow uses:

* A reference linkage map
* A reference VCF file 
* A single chromosome from a reference genome
* Diploid bi-parental F1 population
* Genomic positions for markers order
