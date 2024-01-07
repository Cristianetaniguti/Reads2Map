# 1.5.1

* Make EmpiricalMaps for polyploids compatible with Reads2MapApp v0.0.2
* Polyploid map quality diagnostic by resampling
* Adapt regenotyping memory usage
* Make BarcodeFaker task work for input files that do not finish with .fasta or .fq
* Adapt BarcodeFaker required memory and time
* Adapt STACKs required memory

# 1.5.0

* Adapt tassel and stacks tasks also for polyploids
* Update tests
* Add MAPpoly new functions framework_map and update_framework_map
* Polyploid analysis output compatible with Reads2MapApp v0.0.1
* Remove LargeList deprecated package as dependency

# 1.4.3

* Update example for pair-end reads inputs

# 1.4.2

* more flexibility to choose the probability to be used in the HMM: 

* new parameters:
- global_errors: array with global errors to be tested
- genoprob_error: boolean defining if software probabilities should be tested
- genoprob_global_errors: array with global errors to be tested together with the software probabilities following: 1 - (1 - global error) x (1 - software error probability)

# 1.4.1

* Use BCFtools norm to left-align indel marker positions identified by GATK, STACKs and freebayes

# 1.4.0

* STACKs included
* support to pair-end reads
* defined defaults
* runtimes adapted to run with Caper
* [new tutorial](https://cristianetaniguti.github.io/Tutorials/Reads2Map/Setup_and_run_Reads2Map_workflows.html)

# 1.3.0

* TASSEL 5.0 included
* releases include the input files template
* new parameter to control maximum ram memory used by GATK and TASSEL

# 1.2.0

* Add MAPpoly to build linkage maps for polyploid species
* Run freebayes parallelizing in nodes according to chromosomes and cores splitting in genomic regions
* Adjust runtimes 
* Add polyploid dataset for tests

# 1.0.0

Initial release

This workflow build linkage maps from genotyping-by-sequencing (GBS) data. The GBS samples are splitted into chunks to be run in different nodes and optimize the analyses. Set the number of samples by chunk in the 'chunk_size' input. Use 'max_cores' to define number of cores to be used in each node. The workflow runs the combinations:

* SNP calling: GATK and Freebayes
* Dosage/genotype calling: updog, polyRAD and SuperMASSA
* Linkage map build software: OneMap 3.0 and GUSMap
* Using genotype probabilities from GATK, Freebayes, updog, polyRAD and SuperMASSA, and a global error rate of 5% and 0.001% in the OneMap HMM.

Resulting in 34 linkage maps.

The workflow include de options to:

* Remove or not the read duplicates 
* Perform the Hard Filtering in GATK results
* Replace the VCF AD format field by counts from BAM files
* Run MCHap software to build haplotypes based on GATK called markers
* Include or not multiallelic (MNP) markers
* Apply filters using VCFtools

This workflow requires:

* Diploid bi-parental F1 population
* Single-end reads
* A reference genome
* Genomic positions for markers order
* Selection of a single chromosome from a reference genome to build the linkage maps



