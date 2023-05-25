[![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
[![Reads2Map](https://circleci.com/gh/Cristianetaniguti/Reads2Map.svg?style=svg)](https://app.circleci.com/pipelines/github/Cristianetaniguti/Reads2Map)


<p align="center">
<br>
<img src="https://github.com/Cristianetaniguti/Reads2Map/assets/7572527/6074320a-0eba-44b9-88e1-b89eda8aad70" width="450"/>
<br>
<p/>

Reads2Map presents a collection of [WDL workflows](https://openwdl.org/)  to build linkage maps from sequencing reads. Each workflow release is described in the [Read2Map releases page](https://github.com/Cristianetaniguti/Reads2Map/releases). 

The main workflows are the `EmpiricalReads2Map.wdl` and the `SimulatedReads2Map.wdl`. The `EmpiricalReads2Map.wdl` is composed by the `EmpiricalSNPCalling.wdl` that performs the SNP calling, and the `EmpiricalMaps.wdl` that performs the genotype calling and map building in empirical reads. The `SimulatedReads2Map.wdl` simulates Illumina reads for RADseq, exome, or WGS data and performs the SNP and genotype calling and genetic map building.

By now, [GATK](https://github.com/broadinstitute/gatk), [Freebayes](https://github.com/ekg/freebayes) are included for SNP calling; [updog](https://github.com/dcgerard/updog), [polyRAD](https://github.com/lvclark/polyRAD), [SuperMASSA](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030906) for dosage calling; and [OneMap](https://github.com/augusto-garcia/onemap), [GUSMap](https://github.com/tpbilton/GUSMap), and [MAPpoly](https://github.com/mmollina/MAPpoly) for linkage map build.

The Reads2Map workflows perform the SNP and genotype/dosage calling for your complete data set. However, it builds the linkage map for only a single chromosome (reference genome is required) for each combination of software and parameters. The produced maps will probably still require improvements, but their characteristics will suggest which combination of SNP and genotype calling software and parameters you should use for your data. Once the pipeline is selected, you can input the respective VCF file in R and build the complete linkage map using OneMap or MAPpoly. Use [OneMap](https://statgen-esalq.github.io/tutorials/onemap/Outcrossing_Populations.html) or [MAPoly](https://rpubs.com/mmollin/tetra_mappoly_vignette) tutorials for guidance on building and improving the linkage map for the complete dataset. 

## How to use

Multiple systems are available to run WDL workflows such as Cromwell, miniWDL, and dxWDL. See further information in the [openwdl documentation](https://github.com/openwdl/wdl#execution-engines). 

In addition, we also suggest two wrappers: [cromwell-cli](https://github.com/lmtani/cromwell-cli) and [Caper](https://github.com/ENCODE-DCC/caper). Here is a tutorial on how to setup these tools and one example running the EmpiricalReads2Map:

* [Setup and run Reads2Map workflows](https://cristianetaniguti.github.io/Tutorials/Reads2Map/Setup_and_run_Reads2Map_workflows.html)

To run a pipeline, first navigate to [Reads2Map releases page](https://github.com/Cristianetaniguti/Reads2Map/releases), search for the pipeline tag you which to run, and download the pipeline’s assets (the WDL workflow, the JSON, and the ZIP with accompanying dependencies).

Check the description of the inputs for the pipelines:

* [EmpiricalReads2Map (EmpiricalSNPCalling and EmpiricalMaps)](https://cristianetaniguti.github.io/Tutorials/Reads2Map/EmpiricalReads.html)

* [SimulatedReads2Map](https://cristianetaniguti.github.io/Tutorials/Reads2Map/simulatedreads.html)

Check how to evaluate the workflows results in Reads2MapApp Shiny:

* [Reads2MapApp](https://github.com/Cristianetaniguti/Reads2MapApp)

Check more information and examples of usage in:

* [Taniguti, C. H., Taniguti, L. M., Amadeu, R. R., Mollinari, M., Da, G., Pereira, S., Riera-Lizarazu, O., Lau, J., Byrne, D., de Siqueira Gesteira, G., De, T., Oliveira, P., Ferreira, G. C., &#38; Franco Garcia, A. A.  Developing best practices for genotyping-by-sequencing analysis using linkage maps as benchmarks. BioRxiv. https://doi.org/10.1101/2022.11.24.517847](https://www.biorxiv.org/content/10.1101/2022.11.24.517847v2)

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
- [MCHap](https://github.com/PlantandFoodResearch/MCHap) in [cristaniguti/mchap:0.7.0](https://hub.docker.com/repository/docker/cristaniguti/mchap): Polyploid micro-haplotype assembly using Markov chain Monte Carlo simulation.

### R packages

- [OneMap](https://github.com/augusto-garcia/onemap) in [cristaniguti/reads2map:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/reads2map): Is a software for constructing genetic maps in experimental crosses: full-sib, RILs, F2 and backcrosses;
- [Reads2MapTools](https://github.com/Cristianetaniguti/Reads2MapTools) in [cristaniguti/reads2map:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/reads2map): Support package to perform mapping populations simulations and genotyping for OneMap genetic map building
- [GUSMap](https://github.com/tpbilton/GUSMap): Genotyping Uncertainty with Sequencing data and linkage MAPping
- [updog](https://github.com/dcgerard/updog) in [cristaniguti/reads2map:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/reads2map): Flexible Genotyping of Polyploids using Next Generation Sequencing Data
- [polyRAD](https://github.com/lvclark/polyRAD) in [cristaniguti/reads2map:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/reads2map): Genotype Calling with Uncertainty from Sequencing Data in Polyploids
- [Reads2MapApp](https://github.com/Cristianetaniguti/Reads2MapApp) in [cristaniguti/reads2mapApp:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/reads2map): Shiny app to evaluate Reads2Map workflows results
- [simuscopR](https://github.com/Cristianetaniguti/simuscopR) in [cristaniguti/reads2map:0.0.1](https://hub.docker.com/repository/docker/cristaniguti/reads2map): Wrap-up R package for SimusCop simulations
- [MAPpoly](https://github.com/mmollina/MAPpoly) in [cristaniguti/reads2map:0.0.5](https://hub.docker.com/repository/docker/cristaniguti/reads2map): Build linkage maps for autopolyploid species

### Funding

This work was partially supported by the National Council for Scientific and Technological Development (CNPq - 313269/2021-1); by USDA, National Institute of Food and Agriculture (NIFA), Specialty Crop Research Initiative (SCRI) project “Tools for Genomics Assisted Breeding in Polyploids: Development of a Community Resource” (Award No. 2020-51181-32156); and by the Bill and Melinda Gates Foundation (OPP1213329) project SweetGAINS.