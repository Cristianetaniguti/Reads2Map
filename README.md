## Building linkage maps with onemap_ht

Here you will find workflows written in WDL to perfom linkage map building from simulated and empirical data. The workflow main.wdl provides linkage map simulation using markers from RAD-seq methodology.The workflow empirical.wdl receives fastq files and perform all procedure until the genetic maps. The main goal of the workflow is to measuare the impact of different error probabilities coming from SNP calling methodologies in the linkage map building in OneMap.

## Quickstart

This workflow requires docker hub images. First of all, install [docker](https://docs.docker.com/install/) and [cromwell](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/).

### Run the simulations

* Adapt the path of the inputs in `/main.inputs.json`

*number_of_families* : an integer defining the number of families with `popsize` individuals to be simulated

*references*
- ref_fasta: chromosome sequence in fasta format (only one chromosome at a time)
- ref_fasta_index: index made by samtools faidx
- ref_dict: index made by picard dict
- ref_sa: index made by bwa index
- ref_amb: index made by bwa index
- ref_bwt: index made by bwa index
- ref_ann: index made by bwa index
- ref_pac: index made by bwa index

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

As example, in the folder `data/toy_sample` is available a subset of populus chromosome 10 (find the entire genome [here](https://phytozome.jgi.doe.gov/pz/portal.html#!bulk?org=Org_Ptrichocarpa)), its index files and the `doses` files needed by the workflow. The path to this folder must be defined in `main.inputs.json`.

```
# Execute the workflow
java -jar cromwell.jar run -i main.inputs.json main.wdl
```

**Warning**: See section `WLD Configurations` to choose the better available option for you or create a personalized one.

### Run with empirical data

* Adapt the path of the inputs in `empirical.json`

*empirical.references*
- ref_fasta: chromosome sequence in fasta format (only one chromosome at a time) 
- ref_fasta_index: index made by samtools faidx
- ref_dict: index made by picard dict
- ref_sa: index made by bwa index
- ref_amb: index made by bwa index
- ref_bwt: index made by bwa index
- ref_ann: index made by bwa index
- ref_pac: index made by bwa index

*empirical.dataset*
- samples_info: tsv file with first column with path to fastq file, second column with sample names and third column with sample names and lane specifications.

- name: specie name

#### Test dataset for empirical analysis

You can download black cottonwood genome assembly (FASTA) and RADseq reads from NCBI for testing.

```
# Download a subset of RADseq data from populus study using the docker image "cristaniguti/sratoolkit" with the SRA toolkit 

for i in SRR6249795 SRR6249808 SRR6249768 SRR6249769 SRR6249770 SRR6249771 SRR6249772 SRR6249773 SRR6249774 SRR6249775 SRR6249776 SRR6249778 SRR6249779 SRR6249780 SRR6249781 SRR6249782 SRR6249783 SRR6249784 SRR6249785 SRR6249786 SRR6249787 SRR6249788; do
    docker run -v $(pwd):/opt/ cristaniguti/sratoolkit ./fasterq-dump $i -O /opt/
    head -n 1200000 $i.fastq > $i.sub.fastq # Just a subset of the reads
done

# Here there are enought data to test the pipeline but not for having a good resolution genetic map. It contains the two parents and 20 progeny individuals. The original study have eight replicates for each parent and 122 progenies.
# Access all samples in the BioProject	PRJNA395596

# Download the reference genome in https://phytozome.jgi.doe.gov/pz/portal.html#!bulk?org=Org_Ptrichocarpa 
# Here we only use the first 2500 lines of chromosome 10

```

samples_info_sub

```
data/populus_sub/SRR6249768.sub.fastq   II_3_39 II_3_39.Lib2_B08_TAGGCGC
data/populus_sub/SRR6249769.sub.fastq   II_3_40 II_3_40.Lib2_B09_GGAACTG
data/populus_sub/SRR6249770.sub.fastq   II_3_29 II_3_29.Lib1_A01_AACAATG
data/populus_sub/SRR6249771.sub.fastq   II_3_31 II_3_31.Lib2_B07_ATCTGTT
data/populus_sub/SRR6249772.sub.fastq   II_3_27 II_3_27.Lib2_B05_GGAGACT
data/populus_sub/SRR6249773.sub.fastq   II_3_28 II_3_28.Lib2_B06_CCTCTAG
data/populus_sub/SRR6249774.sub.fastq   II_3_25 II_3_25.Lib2_B03_ACTTCTG
data/populus_sub/SRR6249775.sub.fastq   II_3_26 II_3_26.Lib2_B04_TTGAGGC
data/populus_sub/SRR6249776.sub.fastq   II_3_43 II_3_43.Lib2_B11_TTCTGAG
data/populus_sub/SRR6249778.sub.fastq   I_3_38  I_3_38.Lib1_B10_CCTCAGC
data/populus_sub/SRR6249779.sub.fastq   I_3_34  I_3_34.Lib1_B08_TAGGCGC
data/populus_sub/SRR6249780.sub.fastq   I_3_53  I_3_53.Lib1_C07_TTCTAGT
data/populus_sub/SRR6249781.sub.fastq   I_3_50  I_3_50.Lib1_C06_CCTAGAT
data/populus_sub/SRR6249782.sub.fastq   I_3_45  I_3_45.Lib1_C03_ATGCCGG
data/populus_sub/SRR6249783.sub.fastq   I_3_42  I_3_42.Lib1_C02_CCTGACT
data/populus_sub/SRR6249784.sub.fastq   I_3_59  I_3_59.Lib1_C12_ACGTTGT
data/populus_sub/SRR6249785.sub.fastq   I_3_58  I_3_58.Lib1_C11_TTCCACG
data/populus_sub/SRR6249786.sub.fastq   I_3_56  I_3_56.Lib1_C10_CCTGCAC
data/populus_sub/SRR6249787.sub.fastq   I_3_55  I_3_55.Lib1_C09_AGAAGTC
data/populus_sub/SRR6249788.sub.fastq   I_3_66  I_3_66.Lib1_D06_GCCAACT
data/populus_sub/SRR6249795.sub.fastq   PT_F    PT_F.Lib1_E09_TGAACAT
data/populus_sub/SRR6249808.sub.fastq   PT_M    PT_M.Lib2_E06_CGATGCG
```

```
# Execute the workflow
java -jar cromwell.jar run -i empirical.json empirical.wdl
```

You can also download the full data set running the script "data/populus/download_SRRs.sh" and the "data/populus/sample_info" file.


**Warning**: See section `WLD Configurations` to choose the better available option for you or create a personalized one.

### WDL configurations

There are three possible configurations available in `.configurations` directory:

* cromwell_cache.conf: to store cache in a mysql database

```
# Start a mysql instance on port 3307
# - required just to reuse already processed tasks.
docker run -d -v banco_cromwell:/var/lib/mysql --rm --name mysql-cromwell -p 3307:3306 -e MYSQL_ROOT_PASSWORD=1234 -e MYSQL_DATABASE=cromwell mysql:5.7

# Execute the workflow
java -jar -Dconfig.file=.configurations/cromwell_cache.conf -jar cromwell.jar run -i main.inputs.json main.wdl

```

* cromwell_sing: to run docker images through singularity (useful to run in HPC)

```
# For private repos
export SINGULARITY_DOCKER_USERNAME=fulana
export SINGULARITY_DOCKER_PASSWORD=senhadafulana

# To store singularity temporary files in a different path than home
export SINGULARITY_CACHEDIR=/path/to/cache
export SINGULARITY_TMPDIR=/path/to/cache
export SINGULARITY_LOCALCACHEDIR=/path/to/cache

# Execute the workflow
java -jar -Dconfig.file=.configurations/cromwell_sing.conf -jar cromwell.jar run -i main.inputs.json main.wdl

```

* cromwell_sing_slurm: to run docker images through singularity moderated by slurm system (useful to run in HPC)

```
# Execute the workflow
java -jar -Dconfig.file=.configurations/cromwell_sing_slurm.conf -jar cromwell.jar run -i main.inputs.json main.wdl

```

If you want to run wdl with default configurations simple use:

```
# Execute the workflow
java -jar cromwell.jar run -i main.inputs.json main.wdl
```

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
