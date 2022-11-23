[![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
[![Reads2Map](https://circleci.com/gh/Cristianetaniguti/Reads2Map.svg?style=svg)](https://app.circleci.com/pipelines/github/Cristianetaniguti/Reads2Map)

## Reads2Map 

Reads2Map presents [WDL workflows](https://openwdl.org/) a collection of pipelines to build linkage maps from sequencing reads. Each pipeline release is described in the [Read2Map releases page](https://github.com/Cristianetaniguti/Reads2Map/releases). 

The main workflows are the `EmpiricalSNPCalling.wdl`, the `EmpiricalMaps.wdl`, and the `SimulatedReads.wdl`. `EmpiricalSNPCalling.wdl` performs the SNP calling and `EmpiricalMaps.wdl` performs the genotype calling and map building in empirical reads. The `SimulatedReads.wdl` simulates Illumina reads for RADseq, exome, or WGS data and performs the SNP and genotype calling and genetic map building.

By now, [GATK](https://github.com/broadinstitute/gatk), [Freebayes](https://github.com/ekg/freebayes) are included for SNP calling; [updog](https://github.com/dcgerard/updog), [polyRAD](https://github.com/lvclark/polyRAD), [SuperMASSA](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030906) for dosage calling; and [OneMap](https://github.com/augusto-garcia/onemap), and [GUSMap](https://github.com/tpbilton/GUSMap) for linkage map build.

![math_meth2](https://user-images.githubusercontent.com/7572527/203172239-e4d2d857-84e2-48c5-bb88-01052a287004.png)

## How to use

Multiple systems are available to run WDL workflows such as Cromwell, miniWDL, and dxWDL. See further information in the [openwdl documentation](https://github.com/openwdl/wdl#execution-engines).

To run a pipeline, first navigate to [Reads2Map releases page](https://github.com/Cristianetaniguti/Reads2Map/releases), search for the pipeline tag you which to run, and download the pipeline’s assets (the WDL workflow, the JSON, and the ZIP with accompanying dependencies).

## Documentation

Check the description of the inputs for the pipelines:

* [EmpiricalReads2Map (EmpiricalSNPCalling and EmpiricalMaps)](https://cristianetaniguti.github.io/Tutorials/Read2Map/EmpiricalReads2Map.html)
* [SimulatedReads](https://cristianetaniguti.github.io/Tutorials/Reads2Map/SimulatedReads.html)

Check how to evaluate the workflows results in Reads2MapApp Shiny:

* [Reads2MapApp](https://github.com/Cristianetaniguti/Reads2MapApp)

Once you selected the best pipeline using a subset of your data, you can build a complete high-density linkage map:

* [Quick Guide to build High-Density Linkage Maps](https://cristianetaniguti.github.io/Tutorials/onemap/High_density_maps.html)

Check more information and example of usage in:

* [Paper in preparation]()

## Third-party software and images

- [BWA](https://github.com/lh3/bwa) in [us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z](https://console.cloud.google.com/gcr/images/broad-gotc-prod/US/genomes-in-the-cloud): Used to align simulated reads to reference;
- [cutadapt](https://github.com/marcelm/cutadapt) in [cristaniguti/ pirs-ddrad-cutadapt:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/pirs-ddrad-cutadapt): Trim simulated reads;
- [ddRADseqTools](https://github.com/GGFHF/ddRADseqTools) in [cristaniguti/ pirs-ddrad-cutadapt:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/pirs-ddrad-cutadapt): Set of applications useful to in silico design and testing of double digest RADseq (ddRADseq) experiments;
- [Freebayes](https://github.com/ekg/freebayes) in [Cristaniguti/freebayes:0.0.1](): Variant call step;
- [GATK](https://github.com/broadinstitute/gatk) in [us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z](https://console.cloud.google.com/gcr/images/broad-gotc-prod/US/genomes-in-the-cloud): Variant call step using Haplotype Caller, GenomicsDBImport and GenotypeGVCFs;
- [PedigreeSim](https://github.com/PBR/pedigreeSim?files=1) in [cristaniguti/reads2map:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/reads2map): Simulates progeny genotypes from parents genotypes for different types of populations;
- [picard](https://github.com/broadinstitute/picard) in [us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z](https://console.cloud.google.com/gcr/images/broad-gotc-prod/US/genomes-in-the-cloud): Process alignment files;
- [pirs](https://github.com/galaxy001/pirs) in [cristaniguti/ pirs-ddrad-cutadapt:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/pirs-ddrad-cutadapt): To generate simulates paired-end reads from a reference genome;
- [samtools](https://github.com/samtools/samtools) in [us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z](https://console.cloud.google.com/gcr/images/broad-gotc-prod/US/genomes-in-the-cloud): Process alignment files;
- [SimuSCoP](https://github.com/qasimyu/simuscop) in [cristaniguti/simuscopr:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/simuscopr): Exome and WGS Illumina reads simulations;
- [RADinitio](http://catchenlab.life.illinois.edu/radinitio/) in [	cristaniguti/radinitio:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/radinitio): RADseq Illumina reads simulation;
- [SuperMASSA](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030906) in [cristaniguti/reads2map:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/reads2map): Efficient Exact Maximum a Posteriori Computation for Bayesian SNP Genotyping in Polyploids;
- [bcftools](https://github.com/samtools/bcftools) in [lifebitai/bcftools:1.10.2](https://hub.docker.com/r/lifebitai/bcftools): utilities for variant calling and manipulating VCFs and BCFs;
- [vcftools](http://vcftools.sourceforge.net/) in [cristaniguti/split_markers:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/split_markers): program package designed for working with VCF files.

### R packages

- [OneMap](https://github.com/augusto-garcia/onemap) in [cristaniguti/reads2map:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/reads2map): Is a software for constructing genetic maps in experimental crosses: full-sib, RILs, F2 and backcrosses;
- [Reads2MapTools](https://github.com/Cristianetaniguti/Reads2MapTools) in [cristaniguti/reads2map:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/reads2map): Support package to perform mapping populations simulations and genotyping for OneMap genetic map building
- [GUSMap](https://github.com/tpbilton/GUSMap): Genotyping Uncertainty with Sequencing data and linkage MAPping
- [updog](https://github.com/dcgerard/updog) in [cristaniguti/reads2map:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/reads2map): Flexible Genotyping of Polyploids using Next Generation Sequencing Data
- [polyRAD](https://github.com/lvclark/polyRAD) in [cristaniguti/reads2map:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/reads2map): Genotype Calling with Uncertainty from Sequencing Data in Polyploids
- [Reads2MapApp](https://github.com/Cristianetaniguti/Reads2MapApp) in [cristaniguti/reads2mapApp:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/reads2map): Shiny app to evaluate Reads2Map workflows results
- [simuscopR](https://github.com/Cristianetaniguti/simuscopR) in [cristaniguti/reads2map:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/reads2map): Wrap-up R package for SimusCop simulations.