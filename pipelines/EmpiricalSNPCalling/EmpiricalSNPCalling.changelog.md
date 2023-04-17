# 1.3.0

* TASSEL 5.0 included
* releases include the input files template
* new parameter to control maximum ram memory used by GATK and TASSEL

# 1.2.1

* Adjustments in runtime

# 1.2.0

* Run freebayes parallelizing in nodes according to chromosomes and cores splitting in genomic regions
* Adjust runtimes 
* Add polyploid dataset for tests

# 1.1.0

2022-12-05

* Define default values
* Include optional values to skip some of the sub-workflows
* More tests and example data sets included

# 1.0.0

Initial release

This workflow performs the alignment of FASTQ to a reference genome, SNP calling with GATK tools (HaplotypeCaller, GenomicsDBImport, and GenotypeGVCFs) and Freebayes. The samples are splitted into chunks to be run in different nodes and optimize the analyses. Set the number of samples by chunk in the 'chunk_size' input. Use 'max_cores' to define number of cores to be used in each node.

The workflow also include de options to:

* Remove of not the read duplicates 
* Perform the Hard Filtering in GATK results
* Replace the VCF AD format field by counts from BAM files
* Run MCHap software to build haplotypes based on GATK called markers

This workflow requires:

* Diploid or polyploid specie
* Single-end reads
* A reference genome