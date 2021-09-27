## Reads2Map workflows

Reads2Map workflows offers tools to build linkage maps in diploid outcrossing species from sequencing reads. It compares performances of the selected software: GATK, freebayes, updog, polyRAD, supermassa and test their influences in building genetic maps with OneMap and GUSMap. The main workflows are the `SimulatedReads.wdl` and the `EmpiricalReads.wdl`. The `SimulatedReads.wdl` simulates Illumina reads for RADseq, exome or WGS data and performs the SNP and genotype calling and genetic map building with selected softwares. `EmpiricalReads.wdl` performs the same analysis, but from empirical read sequences.

## Quickstart

These workflows requires docker hub images. First of all, download [Cromwell](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/), install its requirements and [Docker](https://docs.docker.com/install/) or [Singularity](https://sylabs.io/guides/2.6/user-guide/index.html).

### Run SimulatedReads2Map workflow

* Adapt the path of the inputs in `inputs/SimulatedReads2Map.input.json`

*number_of_families* : an integer defining the number of families with `popsize` individuals to be simulated
*global_seed*: This seed is used to generate the families seeds
*max_cores*: Maximum number of computer cores to be used
*filters*: filters in to be applied by VCFtools in the VCF file after SNP calling

*family*: 
- seed: seed to reproduce the analysis after - warning: some steps are still random, as the reads simulation
- popsize: number of individuals at the progenie population
- ploidy: the ploidy of the specie, by now only diploid (2) species are supported
- cross: cross type: by now, only "F1" option is available
- doses: if you do not have an VCF file with variants to be simulated, you can define here the percentage of markers with doses 0, 1 and 2 (when cross is F1)
- cmBymb: if you do not have a reference linkage map, you can simulate using a general recombination rate according with other genetic maps of the specie

*references*
- ref_fasta: chromosome sequence in fasta format (only one chromosome at a time, and no N are allowed)
- ref_fasta_index: index made by samtools faidx
- ref_dict: index made by picard dict
- ref_sa: index made by bwa index
- ref_amb: index made by bwa index
- ref_bwt: index made by bwa index
- ref_ann: index made by bwa index
- ref_pac: index made by bwa index

You can use the docker images to built the indexes files for the reference genome. See an example for `Chr1.2.2M.fa` file presented in the `data/toy_genome/` directory of this repository: 

```
docker run -v $(pwd):/data/ us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z samtools faidx /data/toy_genome/Chr1.2.2M.fa
docker run -v $(pwd):/data/ us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z /usr/gitc/./bwa index /data/toy_genome/Chr1.2.2M.fa
docker run -v $(pwd):/data/ us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z java -jar /usr/gitc/picard.jar CreateSequenceDictionary R=/data/toy_genome/Chr1.2.2M.fa O=/data/toy_genome/Chr1.2.2M.fa
```

or with singularity:

```
singularity run --bind $(pwd):/data/ us.gcr.io_broad-gotc-prod_genomes-in-the-cloud_2.5.7-2021-06-09_16-47-48Z.sif samtools faidx /data/Chr10.populus.fa
singularity run --bind $(pwd):/data/ us.gcr.io_broad-gotc-prod_genomes-in-the-cloud_2.5.7-2021-06-09_16-47-48Z.sif /usr/gitc/./bwa index /data/Chr10.populus.fa
singularity run --bind $(pwd):/data/  us.gcr.io_broad-gotc-prod_genomes-in-the-cloud_2.5.7-2021-06-09_16-47-48Z.sif java -jar /usr/gitc/picard.jar CreateSequenceDictionary R=/data/Chr10.populus.fa O=/data/Chr10.populus.dict
```

The `Chr1.2.2M.fa` contains a subset of [*Populus trichocarpa*](https://phytozome-next.jgi.doe.gov/info/Ptrichocarpa_v4_1) genome.

*sequencing*:
- library_type: the options RADseq, WGS and Exome are available.
- multiallelics: Define with "TRUE" or "FALSE", if the analysis should try to include multiallelic markers in the linkage maps.
- emp_vcf: reference VCF file with the variants to be simulated.
- emp_bam: reference BAM file. It will be used to define the reads profile in WGS and Exome simulation.
- ref_map: reference linkage map, it is a text file with two columns, one named "cM" with values for centimorgan position of markers and other named "bp" with the respective base pair position of each marker. The markers in your reference map do not need to be the same of the VCF file. Using splines, this map is used to train a model to define the position in certimorgan of the simulated variants in the genome.  
- enzyme1: If RADseq, enzyme used the reduce genome representation.
- enzyme2: If RADseq, second enzyme used the reduce genome representation.
- vcf_parent1: parent 1 ID in the reference VCF.
- vcf_parent2: parent 2 ID in the reference VCF.
- chromosome: chromossome ID to be simulated.
- pcr_cycles: If RADseq, the number of PCR cicles used in the library preparation (default: 9).
- insert_size: If RADseq, define the insert size in bp (default: 350).
- read_length: If RADseq, define the read length in bp (default: 150).
- depth: sequencing depth (default: 20).
- insert_size_dev: If RADseq, define the insert size standard deviation in bp (default: 35).


#### Test dataset for simulations

As example, in the directory `data/toy_simulations` you will find input files required to simulate reads and maps based on a subset of *Populus trichocarpa* chromosome 10. These files are: 1) `ref.variants.noindel.recode.vcf` a reference VCF file only with SNPs (indels are not supported by now); 2) and a reference linkage map `ref.map.csv`. The path to the files must be defined in `inputs/SimulatedReads.inputs.json`.

Run the workflow:

```
Execute the workflow
java -jar cromwell.jar run -i SimulatedReads.inputs.json SimulatedReads.wdl
```

**Warning**: See section [Configurations](https://cristianetaniguti.github.io/Tutorials/onemap_workflows/docs/configurations.html) to choose the better available option for you or create a personalized one checking the [cromwell settings for configurations](https://cromwell.readthedocs.io/en/stable/Configuring/).

### Run EmpiricalReads2Map workflow

The EmpiricalReads2Map workflow requires demultiplexed and cleaned FASTQ files. We made available a suggestion for preprocessing reads in `PreprocessingReads.wdl`.

* Adapt the path of the inputs in `inputs/EmpiricalReads2Map.inputs.json`

*empirical.references*
- ref_fasta: chromosome sequence in fasta format (only one chromosome at a time)
- ref_fasta_index: index made by samtools faidx
- ref_dict: index made by picard dict
- ref_sa: index made by bwa index
- ref_amb: index made by bwa index
- ref_bwt: index made by bwa index
- ref_ann: index made by bwa index
- ref_pac: index made by bwa index

You can use docker images to create these indexes, see in [`Run SimulatedReads2Map workflow`](#Run-SimulatedReads2Map-workflow).

*empirical.dataset*
- samples_info: tsv file with first column with path to fastq file, second column with sample names and third column with sample names and lane specifications.

- name: specie name

#### Test dataset for empirical analysis

You can download black cottonwood genome assembly (FASTA) and [RADseq reads from the BioProject	PRJNA395596](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA395596) for testing.

```
# Download a subset of RADseq data from populus study using the docker image "cyverseuk/fastq-dump:latest" 

for i in SRR6249808 SRR6249795 SRR6249788 SRR6249787 SRR6249786 SRR6249785; do
    echo docker run -v $(pwd):/opt/ --rm --entrypoint /usr/bin/fastq-dump -w /opt/ cyverseuk/fastq-dump:latest -X 300000 --gzip $i
done

# Here there are enought data to test the pipeline but not for having a good resolution genetic map. It contains the two parents and 4 progeny individuals. The original study have eight replicates for each parent and 122 progenies.
# Access all samples in the BioProject	PRJNA395596
```

samples_info_sub

```
data/populus_sub/SRR6249785.sub.fastq   I_3_58  I_3_58.Lib1_C11_TTCCACG
data/populus_sub/SRR6249786.sub.fastq   I_3_56  I_3_56.Lib1_C10_CCTGCAC
data/populus_sub/SRR6249787.sub.fastq   I_3_55  I_3_55.Lib1_C09_AGAAGTC
data/populus_sub/SRR6249788.sub.fastq   I_3_66  I_3_66.Lib1_D06_GCCAACT
data/populus_sub/SRR6249795.sub.fastq   PT_F    PT_F.Lib1_E09_TGAACAT
data/populus_sub/SRR6249808.sub.fastq   PT_M    PT_M.Lib2_E06_CGATGCG
```

**Warning**: This analysis demand high computer capacity to run. You will be able to run the example dataset in a computer with 4G of RAM memory, but we suggest to set personalized configurations according to your system. Check some example of configurations in `.configurations` directory. Here is an example of how use them:

BATCH file named `slurm_main.sh`:

```
#!/bin/bash

#SBATCH --export=NONE
#SBATCH -J cromwell_Reads2Map
#SBATCH --nodes=1                    
#SBATCH --mem=1G
#SBATCH --time=01:30:00
#SBATCH -o /home/user/Reads2Map.log
#SBATCH -e /home/user/Reads2Map.err

# Maybe it will be required to import java and singularity modules here. Check the especifications of the HPC.
#module load singularity
#module load java

java -jar -Dconfig.file=/home/user/Reads2Map/.configurations/cromwell_slurm_sing.conf \
     -jar /home/user/Reads2Map/cromwell-65.jar \
     run /home/user/Reads2Map/tasks/snpcalling_emp.wdl \
     -i /home/user/Reads2Map/inputs/toy_sample_emp.inputs.json 
```

Run:

```
sbatch slurm_main.sh
```

Or, for store the metadata and cache in a mysql database (see also [cromwell-cli](https://github.com/lmtani/cromwell-cli) to easy access to metadata):

```
# Open mySQL cointainer
docker run -d -v banco_cromwell:/var/lib/mysql --rm --name mysql-cromwell -p 3307:3306 -e MYSQL_ROOT_PASSWORD=1234 -e MYSQL_DATABASE=cromwell mysql:5.7
# Execute the workflow
java -jar -Dconfig.file=.configurations/cromwell_cache.conf -jar cromwell.jar run -i EmpiricalReads2Map.inputs.json EmpiricalReads2Map.wdl
```

If you want to test the workflow with the full example data, you can download it running "data/populus/download_SRRs.sh" and specify the file "data/populus/sample_info" in the workflow inpu file. In this case and to run full empirical datasets in general, you will need a computer with higher capacity. Check cloud services or if you have access to a High-Performance Computing (HPC). By now, the workflows are configurated for HPC services and we provide some configurations files examples in `.configurations` directory. To adapt it to cloud services it is required to change the `runtime` session of the workflow tasks.

## Visualize Reads2Map workflows output in Reads2MapApp

You can search for all procedure intermediary files in the `cromwell-executions` directory generated by the workflow. The `log` file will specify the workflow id and path for each executed task. The final output is a compressed file called `EmpiricalReads_results.tar.gz` or `SimulatedReads_results.tar.gz`. These files contains tables for an overview of the entire procedure. They are inputs for the Reads2MapApp, an shiny app that provides a graphical view of the results (as presented below). Check the [Reads2MapApp repository](https://github.com/Cristianetaniguti/Reads2MapApp) for further information.

<img width="637" alt="entry" src="https://user-images.githubusercontent.com/7572527/134988930-65e14700-cd13-4b5a-aea8-d27a10402988.PNG">

## Documentation

Here are some tutorials for further details about how to use the workflows:

* [Introduction](https://cristianetaniguti.github.io/Tutorials/onemap_workflows/docs/introduction.html)
* [Running SimulatedReads2Map workflow](https://cristianetaniguti.github.io/Tutorials/onemap_workflows/docs/SimulatedReads2Map.html)
* [Running EmpiricalReads2Map workflow](https://cristianetaniguti.github.io/Tutorials/onemap_workflows/docs/EmpiricalReads2Map.html)
* [High density maps](https://cristianetaniguti.github.io/Tutorials/onemap_workflows/docs/High_density_maps.html)

You can also have more details about the workflows and how they can be applied:

* [Papers]()

## Third party softwares

- [BWA](https://github.com/lh3/bwa): Used to align simulated reads to reference.
- [cleanFastq](https://github.com/davidvi/cleanFastq): Fix broken fastq files with corrupted reads that crash the aligner.
- [cutadapt](https://github.com/marcelm/cutadapt): Trim simulated reads.
- [ddRADseqTools](https://github.com/GGFHF/ddRADseqTools): Set of applications useful to in silico design and testing of double digest RADseq (ddRADseq) experiments.
- [freebayes](https://github.com/ekg/freebayes): Variant call step.
- [GATK](https://github.com/broadinstitute/gatk): Variant call step using Haplotype Caller, GenomicsDBImport and GenotypeGVCFs.
- [OneMap](https://github.com/augusto-garcia/onemap): Is a software for constructing genetic maps in experimental crosses: full-sib, RILs, F2 and backcrosses.
- [GUSMap]()
- [PedigreeSim](https://github.com/PBR/pedigreeSim?files=1): Simulates progeny genotypes from parents genotypes for different types of populations
- [picard](https://github.com/broadinstitute/picard): Process alignment files.
- [pirs](https://github.com/galaxy001/pirs): To generate simulates paired-end reads from a reference genome.
- [samtools](https://github.com/samtools/samtools): Process alignment files.
- [vcf2diploid](https://github.com/abyzovlab/vcf2diploid): Include the variants in a reference genome according with a VCF file.
- [SimuSCoP](https://github.com/qasimyu/simuscop): Exome and WGS Illumina reads simulations.
- [RADinitio](http://catchenlab.life.illinois.edu/radinitio/): RADseq Illumina reads simulation.

