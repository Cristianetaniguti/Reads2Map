# 1.4.0

* STACKs included
* support to pair-end reads
* defined defaults

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



