# 1.0.0

Initial release

This workflow receives as input VCF files from EmpiricalSNPCalling workflow and result in 34 linkage maps for a single chromosome running the combinations:

* SNP calling: GATK and Freebayes
* Dosage/genotype calling: updog, polyRAD and SuperMASSA
* Linkage map build software: OneMap 3.0 and GUSMap
* Using genotype probabilities from GATK, Freebayes, updog, polyRAD and SuperMASSA, and a global error rate of 5% and 0.001% in the OneMap HMM.

It also has the options to:

* Include or not multiallelic (MNP) markers
* Apply filters using VCFtools
