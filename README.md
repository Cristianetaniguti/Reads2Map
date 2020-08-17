## Building linkage maps with onemap_ht

OneMap workflows offers tools to compare performances of the selected software: GATK, freebayes, updog, polyRAD, supermassa and test their influences in building genetic maps with OneMap and GUSMap. The main workflows are the `SimulatedReads.wdl` and the `EmpiricalReads.wdl`. The `SimulatedReads.wdl` simulates Illumina reads for RADseq data and performs the SNP and genotype calling and genetic map building with selected softwares.

## Quickstart

This workflow requires docker hub images. First of all, download [cromwell](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/) and install its requirements and [docker](https://docs.docker.com/install/).

### Run SimulatedReads workflow

* Adapt the path of the inputs in `/SimulatedReads.input.json`

*number_of_families* : an integer defining the number of families with `popsize` individuals to be simulated

*references*
- ref_fasta: chromosome sequence in fasta format (only one chromosome at a time, and no N are allowed)
- ref_fasta_index: index made by samtools faidx
- ref_dict: index made by picard dict
- ref_sa: index made by bwa index
- ref_amb: index made by bwa index
- ref_bwt: index made by bwa index
- ref_ann: index made by bwa index
- ref_pac: index made by bwa index

You can use the docker images to built the indexes files for the reference genome. See a example for `Chr10.populus.fa`:

```
docker run -v $(pwd):/opt/ cristaniguti/r-samtools samtools faidx /opt/Chr10.populus.fa
docker run -v $(pwd):/opt/ kfdrc/bwa-picard:latest-dev bwa index /opt/Chr10.populus.fa
docker run -v $(pwd):/opt/ kfdrc/bwa-picard:latest-dev java -jar picard.jar CreateSequenceDictionary R=/opt/Chr10.populus.fa O=/opt/Chr10.populus.dict
```

*family_template*
- cmBymb: recombination rate according with other genetic maps of the specie
- popsize: number of individuals at the progenie population
- enzyme: enzyme cute site
- seed: seed to reproduce the analysis after - warning: some steps are still random, as the reads simulation
- depth: sequencing depth (remember that the default mode produce pair-end reads, than the double of the depth here defined)
- doses: file containing the percentage of markers with doses 0, 1 and 2 (when cross is F1)
- ploidy: the ploidy of the specie, by now only diploid (2) species are supported
- cross: cross type: "F1" or "F2"

#### Test dataset for simulations

As example, in the folder `data/toy_sample` is available a subset of populus chromosome 10 (find the entire genome [here](https://phytozome.jgi.doe.gov/pz/portal.html#!bulk?org=Org_Ptrichocarpa)), its index files and the `doses` files needed by the workflow. The path to this folder must be defined in `SimulatedReads.inputs.json`.

```
# Execute the workflow
java -jar cromwell.jar run -i SimulatedReads.inputs.json SimulatedReads.inputs.wdl
```

**Warning**: See section [Configurations](documentation/configurations.html) to choose the better available option for you or create a personalized one.

### Run EmpiricalReads workflow

* Adapt the path of the inputs in `EmpiricalReads.inputs.json`

*empirical.references*
- ref_fasta: chromosome sequence in fasta format (only one chromosome at a time)
- ref_fasta_index: index made by samtools faidx
- ref_dict: index made by picard dict
- ref_sa: index made by bwa index
- ref_amb: index made by bwa index
- ref_bwt: index made by bwa index
- ref_ann: index made by bwa index
- ref_pac: index made by bwa index

You can use docker images to create this indexes, see in `Run SimulatedReads workflow`.

*empirical.dataset*
- samples_info: tsv file with first column with path to fastq file, second column with sample names and third column with sample names and lane specifications.

- name: specie name

#### Test dataset for empirical analysis

You can download black cottonwood genome assembly (FASTA) and RADseq reads from NCBI for testing.

```
# Download a subset of RADseq data from populus study using the docker image "cristaniguti/sratoolkit" with the SRA toolkit

for i in SRR6249785 SRR6249786 SRR6249787 SRR6249788; do
    docker run -v $(pwd):/opt/ cristaniguti/sratoolkit ./fasterq-dump $i -O /opt/
    head -n 80000 $i.fastq > $i.sub.fastq # Just a subset of the reads
done

# Here there are enought data to test the pipeline but not for having a good resolution genetic map. It contains the two parents and 20 progeny individuals. The original study have eight replicates for each parent and 122 progenies.
# Access all samples in the BioProject	PRJNA395596

# Download the reference genome in https://phytozome.jgi.doe.gov/pz/portal.html#!bulk?org=Org_Ptrichocarpa
# Here we only use the first 2500 lines of chromosome 10

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

**Warning**: This analysis demand more computer capacity to run. Then, we suggest you choose a configuration in [Configurations](documentation/configurations.html) before run, as below.

```
# Open mySQL cointainer
docker run -d -v banco_cromwell:/var/lib/mysql --rm --name mysql-cromwell -p 3307:3306 -e MYSQL_ROOT_PASSWORD=1234 -e MYSQL_DATABASE=cromwell mysql:5.7
# Execute the workflow
java -jar -Dconfig.file=.configurations/cromwell_cache.conf -jar cromwell.jar run -i main.inputs.json main.wdl
```

You can also download the full data set running the script "data/populus/download_SRRs.sh" and run the workflow using "data/populus/sample_info" file.


## Documentation

Here are some tutorials that better explain how to use the workflows:

* [Introduction](https://cristianetaniguti.github.io/onemap_workflows/introduction.html)
* [Running SimulatedReads workflow](https://cristianetaniguti.github.io/onemap_workflows/simulatedreads.html)
* [Running EmpiricalReads workflow](https://cristianetaniguti.github.io/onemap_workflows/empiricalreads.html)
* [Configurations](https://cristianetaniguti.github.io/onemap_workflows/configuration.html)
* [Usage of sub-workflows](https://cristianetaniguti.github.io/onemap_workflows/subworkflows.html)

You can also have more details about the workflows and how they can be applied:

* [Technical note]()
* [Paper errors]()
* [Paper MNPs]()

## Third party softwares

- [bwa](https://github.com/lh3/bwa): Used to align simulated reads to reference.
- [cleanFastq](https://github.com/davidvi/cleanFastq): Fix broken fastq files with corrupted reads that crash the aligner.
- [cutadapt](https://github.com/marcelm/cutadapt): Trim simulated reads.
- [ddRADseqTools](https://github.com/GGFHF/ddRADseqTools): Set of applications useful to in silico design and testing of double digest RADseq (ddRADseq) experiments.
- [freebayes](https://github.com/ekg/freebayes): Variant call step.
- [gatk](https://github.com/broadinstitute/gatk): Variant call step using Haplotype Caller, GenomicsDBImport and GenotypeGVCFs.
- [onemap](https://github.com/augusto-garcia/onemap): Is a software for constructing genetic maps in experimental crosses: full-sib, RILs, F2 and backcrosses.
- [PedigreeSim](https://github.com/PBR/pedigreeSim?files=1): Simulates progeny genotypes from parents genotypes for different types of populations
- [picard](https://github.com/broadinstitute/picard): Process alignment files.
- [pirs](https://github.com/galaxy001/pirs): To generate simulates paired-end reads from a reference genome.
- [samtools](https://github.com/samtools/samtools): Process alignment files.
- [vcf2diploid](https://github.com/abyzovlab/vcf2diploid): Include the variants in a reference genome according with a VCF file.
