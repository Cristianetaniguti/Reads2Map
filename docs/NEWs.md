# onemap 3.0

* split_onemap function to split onemap objects
* onemap_read_vcfR now works for multiallelic markers
* map parallelized according with Batchmap method
* parmap function - alternative function to parallelize HMM
* HMM now uses three possible errors in emission phase: global error, genotype errors, genotype probabilities
* create_probs function makes the convertion of the three types of errors to the emission matrix
* updog_genotype function performs regenotyping with updog software
* polyRAD_genotype function performs regenotyping with polyRAD software
* create_depths_profile plots allele counts and genotypes from vcf, onemap and errors
* runpedsim function makes a interface with PedigreeSim software
* pedsim2raw converts PedigreeSim outputs to onemap raw file
* pedsim2vcf converts PedigreeSim outputs to VCF file simulating counts with updog or negative binomial
* Replace all get() by nothing
* New vignettes "Simulations" and "High Density Maps"
* version compatible with onemap_workflows and onemap_workflows_app shiny app

Add it in README

5. [How to simulate maps](http://critianetaniguti.github.io/onemap/vignettes_highres/Simulations.html)

6. [How to build a high density  maps with markers from high-throughput sequencing](http://critianetaniguti.github.io/onemap/vignettes_highres/High_density_maps.html)