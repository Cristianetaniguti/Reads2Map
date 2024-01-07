# 1.3.1

* Make EmpiricalMaps for polyploids compatible with Reads2MapApp v0.0.2
* Polyploid map quality diagnostic by resampling
* Adapt regenotyping memory usage

# 1.3.0

* Add MAPpoly new functions framework_map and update_framework_map
* Update tests
* Polyploid analysis output compatible with Reads2MapApp v0.0.1
* Remove LargeList deprecated package as dependency

# 1.2.5

* more flexibility to choose the probability to be used in the HMM:

* new parameters:
- global_errors: array with global errors to be tested
- genoprob_error: boolean defining if software probabilities should be tested
- genoprob_global_errors: array with global errors to be tested together with the software probabilities following: 1 - (1 - global error) x (1 - software error probability)

# 1.2.4

* runtimes adapted to run with Caper
* perform the genotype calling with updog, SuperMASSA and polyRAD with complete data set (not only for the selected chromosome)
* [new tutorial](https://cristianetaniguti.github.io/Tutorials/Reads2Map/Setup_and_run_Reads2Map_workflows.html)

# 1.2.3

* Supermassa has smaller probability threshold (bugfix)

# 1.2.2

* Supermassa has smaller probability threshold

# 1.2.1

* Avoid estimating multipoint genetic distances to save time
* Adjustments in runtime

# 1.2.0

* Add MAPpoly to build linkage maps for polyploid species
* Adjust runtimes 
* Add polyploid dataset for tests

# 1.1.0

2022-12-05

* Include optional values to skip some of the sub-workflows
* Accepts any number of VCF files as input
* Dosage calling with updog, polyRAD and SuperMASSA adapted to polyploids
* More tests and example data sets included

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

This workflow uses:

* Diploid bi-parental F1 population
* Genomic positions for markers order
* A single chromosome from a reference genome
