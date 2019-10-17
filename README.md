## Linkage map simulation

This workflow provides linkage map simulation using markers from RAD-seq methodology. The main goal of the workflow is to measuare the impact of different error probabilities coming from SNP calling methodologies in the linkage map building in OneMap.

## Quickstart

This workflow requires docker hub images. First of all, install [docker](https://docs.docker.com/install/) and [cromwell](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/).

### Run for only one family

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
