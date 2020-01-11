version 1.0

import "structs/empiricalS.wdl"

workflow empirical{

 input{
   dataset dataset
   ReferenceFasta references
 }

 Array[Array[File]] inputSamples = read_tsv(dataset.samples_info)

 scatter (samples in inputSamples){
   call RunBwaAlignment{
     input:
       sampleName = samples[2],
       reads1     = samples[0],
       ref        = references.ref_fasta,
       geno_amb   = references.ref_amb,
       geno_ann   = references.ref_ann,
       geno_bwt   = references.ref_bwt,
       geno_pac   = references.ref_pac,
       geno_sa    = references.ref_sa
   }

   call AddAlignmentHeader{
     input:
       sampleName = samples[1],
       libName    = samples[2],
       bam_file   = RunBwaAlignment.bam_file,
       bam_idx    = RunBwaAlignment.bam_idx
   }
 }


 call JointSameSamples{
   input:
     samples_info = dataset.samples_info,
     bam_rg       = AddAlignmentHeader.bam_rg
 }

 Array[String] merged_names = read_lines(JointSameSamples.merged_names)
 Array[Pair[File, String]] bam_files = zip(JointSameSamples.merged_files, merged_names)

 scatter (bams in bam_files){

   call HaplotypeCallerERC {
     input:
       ref        = references.ref_fasta,
       geno_fai   = references.ref_fasta_index,
       sampleName = bams.right,
       bam_rg     = bams.left,
       bam_rg_idx = JointSameSamples.merged_files_idx,
       geno_dict  = references.ref_dict
   }

 }

 call CreateGatkDatabase {
   input:
     path_gatkDatabase = "my_database",
     GVCFs             = HaplotypeCallerERC.GVCF,
     GVCFs_idx         = HaplotypeCallerERC.GVCF_idx
 }

 call GenotypeGVCFs {
   input:
     workspace_tar       = CreateGatkDatabase.workspace_tar,
     output_vcf_filename = dataset.name + "_gatk.vcf",
     ref                 = references.ref_fasta,
     geno_fai            = references.ref_fasta_index,
     geno_dict           = references.ref_dict
 }


 call RunFreebayes {
   input:
     freebayesVCFname = dataset.name + "_freebayes.vcf",
     ref              = references.ref_fasta,
     ref_fai          = references.ref_fasta_index,
     merged_files     = JointSameSamples.merged_files
 }

 call VcftoolsApplyFilters{
   input:
     freebayesVCF = RunFreebayes.freebayesVCF,
     gatkVCF      = GenotypeGVCFs.gatkVCF
 }

 scatter (bams in bam_files) {
   call BamCounts{
     input:
       sampleName     = bams.right,
       bam_file       = bams.left,
       bam_idx        = JointSameSamples.merged_files_idx,
       ref            = references.ref_fasta,
       ref_fai        = references.ref_fasta_index,
       ref_dict       = references.ref_dict,
       gatk_vcf       = VcftoolsApplyFilters.gatkVCF_F,
       freebayes_vcf  = VcftoolsApplyFilters.freebayesVCF_F
   }
 }

 call BamCounts4Onemap {
   input:
     samples_info       = dataset.samples_info,
     freebayes_counts = BamCounts.freebayes_counts,
     gatk_counts      = BamCounts.gatk_counts
 }

 Array[String] methods = ["gatk", "freebayes"]
 Array[File] vcfs = [VcftoolsApplyFilters.gatkVCF_F, VcftoolsApplyFilters.freebayesVCF_F]
 Array[Pair[String, File]] program_and_vcf = zip(methods, vcfs)

 scatter (vcf in program_and_vcf){
   call avalVCFs{
     input:
       methodName = vcf.left,
       vcf_file   = vcf.right,
       freebayes_ref_depth = BamCounts4Onemap.freebayes_ref_bam,
       freebayes_alt_depth = BamCounts4Onemap.freebayes_alt_bam,
       gatk_ref_depth = BamCounts4Onemap.gatk_ref_bam,
       gatk_alt_depth = BamCounts4Onemap.gatk_alt_bam,
       parent1 = dataset.parent1,
       parent2 = dataset.parent2,
       gatk_example_alleles = BamCounts4Onemap.gatk_example_alleles,
       freebayes_example_alleles = BamCounts4Onemap.freebayes_example_alleles
   }
 }

  call JointDatas{
      input:
       depths_df_vcf_rds = avalVCFs.depths_df_vcf_rds,
       depths_updog_vcf_rds = avalVCFs.depths_updog_vcf_rds,
       depths_supermassa_vcf_rds = avalVCFs.depths_supermassa_vcf_rds,
       depths_polyrad_vcf_rds = avalVCFs.depths_polyrad_vcf_rds,
       depths_updog_bam_rds = avalVCFs.depths_updog_bam_rds,
       depths_supermassa_bam_rds = avalVCFs.depths_supermassa_bam_rds,
       depths_polyrad_bam_rds = avalVCFs.depths_polyrad_bam_rds,
       times_vcf_df = avalVCFs.times_vcf_df,
       times_vcf_GQ = avalVCFs.times_vcf_GQ,
       times_vcf_updog = avalVCFs.times_vcf_updog,
       times_vcf_supermassa = avalVCFs.times_vcf_supermassa,
       times_vcf_polyrad = avalVCFs.times_vcf_polyrad,
       times_vcf_gusmap = avalVCFs.times_vcf_gusmap,
       times_bam_gusmap = avalVCFs.times_bam_gusmap,
       times_bam_updog = avalVCFs.times_bam_updog,
       times_bam_supermassa = avalVCFs.times_bam_supermassa,
       times_bam_polyrad = avalVCFs.times_bam_polyrad,
       filters_vcf_dfAndGQ = avalVCFs.filters_vcf_dfAndGQ,
       filters_vcf_polyrad = avalVCFs.filters_vcf_polyrad,
       filters_vcf_supermassa = avalVCFs.filters_vcf_supermassa,
       filters_vcf_updog = avalVCFs.filters_vcf_updog,
       filters_bam_polyrad = avalVCFs.filters_bam_polyrad,
       filters_bam_supermassa = avalVCFs.filters_bam_supermassa,
       filters_bam_updog = avalVCFs.filters_bam_updog,
       RData_vcf_df = avalVCFs.RData_vcf_df,
       RData_vcf_GQ = avalVCFs.RData_vcf_GQ,
       RData_vcf_polyrad = avalVCFs.RData_vcf_polyrad,
       RData_vcf_supermassa = avalVCFs.RData_vcf_supermassa,
       RData_vcf_updog = avalVCFs.RData_vcf_updog,
       RData_bam_polyrad = avalVCFs.RData_bam_polyrad,
       RData_bam_supermassa = avalVCFs.RData_bam_supermassa,
       RData_bam_updog = avalVCFs.RData_bam_updog,
       RData_gusmap = avalVCFs.RData_gusmap,
       RData_bam_gusmap = avalVCFs.RData_bam_gusmap,
       map_vcf_df = avalVCFs.map_vcf_df,
       map_vcf_GQ = avalVCFs.map_vcf_GQ,
       map_vcf_polyrad = avalVCFs.map_vcf_polyrad,
       map_vcf_supermassa = avalVCFs.map_vcf_supermassa,
       map_vcf_updog = avalVCFs.map_vcf_updog,
       map_bam_polyrad = avalVCFs.map_bam_polyrad,
       map_bam_supermassa = avalVCFs.map_bam_supermassa,
       map_bam_updog = avalVCFs.map_bam_updog,
       map_gusmap = avalVCFs.map_gusmap,
       map_bam_gusmap = avalVCFs.map_bam_gusmap
       }

       output{
               File data_depths = JointDatas.data_depths
               File data_filters = JointDatas.data_filters
               File data_map = JointDatas.data_map
               File data_times = JointDatas.data_times
               File data_RDatas = JointDatas.data_RDatas
       }
}

task RunBwaAlignment {

 input {
   String sampleName
   File ref
   File reads1
   File geno_amb
   File geno_ann
   File geno_bwt
   File geno_pac
   File geno_sa
 }

 command <<<
   export PATH=$PATH:/bin
   export PATH=$PATH:/picard.jar

   bwa mem -t 10 ~{ref} ~{reads1} | \
     java -jar /picard.jar SortSam \
       I=/dev/stdin \
       O=~{sampleName}.sorted.bam \
       SORT_ORDER=coordinate \
       CREATE_INDEX=true
   mv ~{sampleName}.sorted.bai ~{sampleName}.sorted.bam.bai
 >>>

 runtime {
   docker: "kfdrc/bwa-picard:latest-dev"
 }

 output {
   File bam_file = "${sampleName}.sorted.bam"
   File bam_idx = "${sampleName}.sorted.bam.bai"
 }
}

# Add info to alignment header
task AddAlignmentHeader {
 input {
   String sampleName
   String libName
   File bam_file
   File bam_idx
 }

 command <<<
   mkdir tmp
   java -jar /gatk/picard.jar AddOrReplaceReadGroups \
     I=~{bam_file} \
     O=~{libName}_rg.bam \
     RGLB=lib-~{libName} \
     RGPL=illumina \
     RGID=FLOWCELL1.LANE1.~{libName} \
     RGSM=~{sampleName} \
     RGPU=FLOWCELL1.LANE1.~{libName} \
     CREATE_INDEX=true \
     TMP_DIR=tmp

   mv ~{libName}_rg.bai ~{libName}_rg.bam.bai
 >>>

 runtime {
   docker: "taniguti/gatk-picard"
 }

 output {
   File bam_rg = "${libName}_rg.bam"
   File bam_rg_index = "${libName}_rg.bam.bai"
 }
}

# Joint same sample bam files
task JointSameSamples{

 input{
   File samples_info
   Array[File] bam_rg
 }

 command <<<

   R --vanilla --no-save <<RSCRIPT

     system("cp ~{sep=" "  bam_rg} .")

     files <- read.table("~{samples_info}", stringsAsFactors = F)

     repet <- names(which(table(files[,2]) > 1))

     if(length(repet) != 0){
       idx <- vector()
       for(i in 1:length(repet)){
         idx <- c(idx,which(files[,2] == repet[i]))
         files1 <- files[which(files[,2] == repet[i]),3]
         files1 <- paste0(files1, "_rg.bam")
         system(paste0("samtools merge ", repet[i], ".merged.bam"," ", paste(files1, collapse = " ") , collapse=" "))
       }
       files2 <- files[-idx,]
     } else {
       files2 <- files
     }

     for(i in 1:dim(files2)[1]){
       system(paste0("mv ", files2[,3][i], "_rg.bam ", files2[,2][i], ".merged.bam "))
     }

     system("ls *merged.bam > merged_names")
     df <- read.table("merged_names")
     for(i in 1:length(df[,1]))
         system(paste0("samtools index ", df[i,1]))
     df.new <- sapply(strsplit(as.character(df[,1]), "[.]"), "[",1)
     write.table(df.new, "merged_names", quote = F, col.names=F, row.names=F)

   RSCRIPT
 >>>

 runtime{
   docker: "cristaniguti/r-samtools"
 }

 output{
   Array[File] merged_files = glob("*.merged.bam")
   Array[File] merged_files_idx = glob("*.merged.bam.bai")
   File merged_names = "merged_names"
 }
}

# GATK to generate gVCF with variants
task HaplotypeCallerERC {
 input {
   File ref
   File geno_fai
   String sampleName
   File bam_rg
   Array[File] bam_rg_idx
   File geno_dict
 }

 command <<<

   cp ~{sep=" " bam_rg_idx} $(dirname ~{bam_rg})

   /gatk/gatk HaplotypeCaller \
     -ERC GVCF \
     -R ~{ref} \
     -I ~{bam_rg} \
     -O ~{sampleName}_rawLikelihoods.g.vcf \
     --max-reads-per-alignment-start 0
 >>>

 runtime {
   docker: "taniguti/gatk-picard"
 }

 output {
   File GVCF = "${sampleName}_rawLikelihoods.g.vcf"
   File GVCF_idx = "${sampleName}_rawLikelihoods.g.vcf.idx"
 }
}

task CreateGatkDatabase {

 input {
   String path_gatkDatabase
   Array[File] GVCFs
   Array[File] GVCFs_idx
 }

 command <<<
   /gatk/gatk GenomicsDBImport \
     --genomicsdb-workspace-path ~{path_gatkDatabase} \
     -L Chr10 \
     -V ~{sep=" -V "  GVCFs}

   tar -cf ~{path_gatkDatabase}.tar ~{path_gatkDatabase}
 >>>

 runtime {
   docker: "taniguti/gatk-picard"
 }

 output {
   File workspace_tar = "${path_gatkDatabase}.tar"
 }
}

# Variant calling on gVCF
task GenotypeGVCFs {

 input {
   File workspace_tar
   String output_vcf_filename
   File ref
   File geno_fai
   File geno_dict
 }

 command <<<
   tar -xf ~{workspace_tar}
   WORKSPACE=$( basename ~{workspace_tar} .tar)

   /gatk/gatk GenotypeGVCFs \
       -R ~{ref} \
       -O ~{output_vcf_filename} \
       -G StandardAnnotation \
       -V gendb://$WORKSPACE
 >>>

 runtime {
   docker: "taniguti/gatk-picard"
 }

 output {
   File gatkVCF = "${output_vcf_filename}"
   File gatkVCF_index = "${output_vcf_filename}.idx"
 }
}


# Variant calling using freebayes
task RunFreebayes {

 input {
   String freebayesVCFname
   File ref
   File ref_fai
   Array[File] merged_files
 }

 command <<<
   freebayes --genotype-qualities -f ~{ref} ~{sep=" "  merged_files} > ~{freebayesVCFname}
 >>>

 runtime {
   docker: "taniguti/freebayes"
 }

 output {
   File freebayesVCF = "${freebayesVCFname}"
 }
}


task  VcftoolsApplyFilters {

 input {
   File gatkVCF
   File freebayesVCF
 }

 command <<<
   vcftools --vcf "~{gatkVCF}" --max-missing 0.75 --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --out gatk
   vcftools --vcf "~{freebayesVCF}" --max-missing 0.75 --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --out freebayes

 >>>
 runtime {
   docker: "taniguti/vcftools"
 }

 output {
   File gatkVCF_F = "gatk.recode.vcf"
   File freebayesVCF_F = "freebayes.recode.vcf"
 }
}


# This task extract the allele depths from bam files
task BamCounts{
 input{
   String sampleName
   File bam_file
   Array[File] bam_idx
   File ref
   File ref_fai
   File ref_dict
   File gatk_vcf
   File freebayes_vcf
 }

 command <<<

   cp ~{sep=" " bam_idx} $(dirname ~{bam_file})

   java -jar /gatk/picard.jar VcfToIntervalList \
     I=~{gatk_vcf} \
     O=gatk.interval.list

   /gatk/gatk CollectAllelicCounts \
     --input ~{bam_file} \
     --reference ~{ref} \
     --intervals gatk.interval.list \
     --output ~{sampleName}_gatk_counts.tsv

  java -jar /gatk/picard.jar VcfToIntervalList \
     I=~{freebayes_vcf} \
     O=freebayes.interval.list

   /gatk/gatk CollectAllelicCounts \
     --input ~{bam_file} \
     --reference ~{ref} \
     --intervals freebayes.interval.list \
     --output ~{sampleName}_freebayes_counts.tsv

 >>>

 runtime{
   docker:"taniguti/gatk-picard"
 }

 output{
   File gatk_counts = "~{sampleName}_gatk_counts.tsv"
   File freebayes_counts = "~{sampleName}_freebayes_counts.tsv"
 }
}

task BamCounts4Onemap{
 input{
   Array[File] freebayes_counts
   Array[File] gatk_counts
   File samples_info
 }

 command <<<
   R --vanilla --no-save <<RSCRIPT
     library(R.utils)
     system("cp ~{sep=" "  freebayes_counts} .")
     system("cp ~{sep=" "  gatk_counts} .")
     info <- read.table("~{samples_info}", stringsAsFactors=F)
     names <- info[,2]
     methods <- c("freebayes", "gatk")

     for(method in methods){

     file.counts <- read.table(paste0(names[1],"_", method,"_counts.tsv"), skip = 3, header=T, stringsAsFactors = F)

     ref_depth_matrix2 <- alt_depth_matrix2  <- matrix(NA, nrow = dim(file.counts)[1], ncol = length(names))

     for(j in 1:length(names)){
       ## From picard tool

       file.counts <- read.table(paste0(names[j],"_", method,"_counts.tsv"), skip = 3, header=T, stringsAsFactors = F)

       ref_depth_matrix2[,j] <- file.counts[,3]
       alt_depth_matrix2[,j] <- file.counts[,4]

       if(j == 1){
         ref_allele <- file.counts[,5]
         alt_allele <- file.counts[,6]
       } else {
         idx.ref <- which(ref_allele == "N")
         idx.alt <- which(alt_allele == "N")
         if(length(idx.ref)!=0){
           ref_allele[idx.ref] <- file.counts[idx.ref,5]
         }
         if(length(idx.alt)!=0){
           alt_allele[idx.alt] <- file.counts[idx.alt,6]
         }
       }

     }

     rownames(ref_depth_matrix2) <- rownames(alt_depth_matrix2) <- paste0(file.counts[,1],"_", file.counts[,2])
     colnames(ref_depth_matrix2) <- colnames(alt_depth_matrix2) <- names

     alleles <- data.frame(file.counts[,1],file.counts[,2], ref_allele, alt_allele)
     write.table(alleles, file = paste0(method,"_example4ref_alt_alleles.txt"), col.names = F, row.names = F)

     write.table(ref_depth_matrix2, file = paste0(method,"_ref_depth_bam.txt"), quote=F, row.names=T, sep="\t", col.names=T)
     write.table(alt_depth_matrix2, file = paste0(method,"_alt_depth_bam.txt"), quote=F, row.names=T, sep="\t", col.names=T)
   }

   RSCRIPT
 >>>

 runtime{
   docker:"taniguti/onemap"
 }

 output{
   File gatk_ref_bam      = "gatk_ref_depth_bam.txt"
   File gatk_alt_bam      = "gatk_alt_depth_bam.txt"
   File freebayes_ref_bam = "freebayes_ref_depth_bam.txt"
   File freebayes_alt_bam = "freebayes_alt_depth_bam.txt"
   File gatk_example_alleles    = "gatk_example4ref_alt_alleles.txt"
   File freebayes_example_alleles    = "freebayes_example4ref_alt_alleles.txt"
 }
}


task avalVCFs{
 input{
   String methodName
   File vcf_file
   File freebayes_ref_depth
   File freebayes_alt_depth
   File gatk_ref_depth
   File gatk_alt_depth
   String parent1
   String parent2
   File gatk_example_alleles
   File freebayes_example_alleles
 }

 command <<<

   R --vanilla --no-save <<RSCRIPT

   system("cp ~{freebayes_ref_depth} .")
   system("cp ~{freebayes_alt_depth} .")
   system("cp ~{gatk_ref_depth} .")
   system("cp ~{gatk_alt_depth} .")
   system("cp ~{gatk_example_alleles} .")
   system("cp ~{freebayes_example_alleles} .")

   source("/opt/scripts/functions_empirical.R")

   vcf_file <- "~{vcf_file}"
   SNPCall <- "~{methodName}"
   cross <- "outcross"
   max.cores <- 3
   parent1 <- "~{parent1}"
   parent2 <- "~{parent2}"

   # PACKAGES
   library(supermassa4onemap)
   library(onemap)
   library(updog)
   library(reshape2)
   library(vcfR)
   library(doParallel)
   library(GUSMap)
   library(ggplot2)

   ## READING FINAL VCF FROM PIPELINE
   vcf <- read.vcfR(vcf_file)
   df <- onemap_read_vcfR(vcfR.object=vcf,
                         cross= cross,
                         parent1= parent1,
                         parent2= parent2)

   # check depths
   p <- create_depths_profile(onemap.obj = df, vcfR.object = vcf, parent1 = parent1, parent2 = parent2, vcf.par = "AD", recovering = FALSE, GTfrom = "vcf", alpha=0.1,
   rds.file = paste0(SNPCall,"_", GenoCall= "df","_",CountsFrom="vcf","_depths.rds"))
   #ggsave(filename = paste0(SNPCall,"_", GenoCall= "df","_",CountsFrom="vcf","_vcf_depths.png"), p)

   ## FILTERS REPORT
   out_name <- paste0(SNPCall, "_filters_vcf_dfAndGQ.txt")
   filters_tab <- create_filters_report(onemap_obj = df, CountsFrom = "vcf", SNPCall=SNPCall, GenoCall="df")
   write_report(filters_tab[[1]], out_name)

   # MAPS REPORT - DF
   create_map_report(input.seq = filters_tab[[2]], CountsFrom = "vcf", SNPCall = SNPCall, GenoCall="df")

   # MAPS REPORT - GQ
   aval.gq <- extract_depth(vcfR.object=vcf,
                           onemap.object=df,
                           vcf.par="GQ",
                           parent1=parent1,
                           parent2=parent2,
                           recovering=FALSE)

   aval.gq <- create_probs(df, genotypes_errors=aval.gq)
   filters_tab <- create_filters_report(aval.gq, CountsFrom = "vcf", SNPCall=SNPCall, GenoCall="df")

   create_map_report(filters_tab[[2]], CountsFrom = "vcf", SNPCall = SNPCall, GenoCall="GQ")

   # OTHER TOOLS
   ## With depths from vcf

   updog.aval <- updog_error(
     vcfR.object=vcf,
     onemap.object=df,
     vcf.par="AD",
     parent1=parent1,
     parent2=parent2,
     recovering=TRUE,
     mean_phred=20,
     cores=max.cores,
     depths=NULL)

   supermassa.aval <- supermassa4onemap::supermassa_error(
     vcfR.object=vcf,
     onemap.object = df,
     vcf.par = "AD",
     parent1 = parent1,
     parent2 = parent2,
     recovering = TRUE,
     mean_phred = 20,
     cores = max.cores,
     depths = NULL)

   polyrad.aval <- polyRAD_error(
     vcf=vcf_file,
     onemap.obj=df,
     parent1=parent1,
     parent2=parent2,
     crosstype=cross)

   metodologies <- list(updog = updog.aval, supermassa = supermassa.aval, polyrad = polyrad.aval)
   for (metodology in names(metodologies)){
     cat(metodology, "\n")
     error_aval <- metodologies[[metodology]]
     ## Filters
     out_name <- paste0(SNPCall, "_filters_vcf_", metodology, ".txt")
     filters_tab <- create_filters_report(error_aval, CountsFrom = "vcf", SNPCall = SNPCall, GenoCall = metodology)
     write_report(filters_tab[[1]], out_name)
     
     # check depths
     p <- create_depths_profile(onemap.obj = error_aval, vcfR.object = vcf, parent1 = parent1, parent2 = parent2, vcf.par = "AD", 
                               recovering = FALSE, GTfrom = "onemap", alpha=0.1,
                               rds.file = paste0(SNPCall,"_", GenoCall= metodology, "_",CountsFrom = "vcf","_depths.rds"))
                       
     #ggsave(filename = paste0(SNPCall,"_", GenoCall= metodology, "_",CountsFrom = "vcf","_onemap_depths.png"), p)
     
     ## Maps
     create_map_report(input.seq = filters_tab[[2]], CountsFrom = "vcf", SNPCall = SNPCall, GenoCall= metodology)
   }

   ## Depths from bam
   depths.alt <- read.table(paste0(SNPCall, "_alt_depth_bam.txt"), header = T)
   depths.ref <- read.table(paste0(SNPCall, "_ref_depth_bam.txt"), header = T)

   depths <- list("ref"=depths.ref, "alt"=depths.alt)

   updog.aval.bam <- updog_error(
     vcfR.object=vcf,
     onemap.object=df,
     vcf.par="AD",
     parent1=parent1,
     parent2=parent2,
     recovering=TRUE,
     mean_phred=20,
     cores=max.cores,
     depths=depths)

   supermassa.aval.bam <- supermassa_error(
     vcfR.object=vcf,
     onemap.object = df,
     vcf.par = "AD",
     parent1 = parent1,
     parent2 = parent2,
     recovering = TRUE,
     mean_phred = 20,
     cores = max.cores,
     depths = depths)

   new.vcf <- make_vcf(vcf_file, depths, SNPCall)
   new.vcfR <- read.vcfR(new.vcf) 

   polyrad.aval.bam <- polyRAD_error(
     vcf=new.vcf,
     onemap.obj=df,
     parent1=parent1,
     parent2=parent2,
     crosstype=cross)

   metodologies <- list(updog = updog.aval.bam, supermassa= supermassa.aval.bam, polyrad=polyrad.aval.bam)
   for (metodology in names(metodologies)){
     cat(metodology, "\n")
     error_aval <- metodologies[[metodology]]
     ## Filters
     out_name <- paste0(SNPCall, "_filters_bam_", metodology, ".txt")
     filters_tab <- create_filters_report(error_aval, CountsFrom = "bam", SNPCall = SNPCall, GenoCall = metodology)
     write_report(filters_tab[[1]], out_name)

     # check depths
     p <- create_depths_profile(onemap.obj = error_aval, vcfR.object = new.vcfR, parent1 = parent1, parent2 = parent2, vcf.par = "AD", 
                               recovering = FALSE, GTfrom = "onemap", alpha=0.1,
                               rds.file = paste0(SNPCall,"_", GenoCall= metodology, "_",CountsFrom = "bam","_depths.rds"))
     #ggsave(filename = paste0(SNPCall,"_", GenoCall= metodology, "_",CountsFrom = "bam","_onemap_depths.png"), p)
     
     ## Maps
     create_map_report(input.seq = filters_tab[[2]], CountsFrom = "bam", SNPCall = SNPCall, GenoCall= metodology)
   }

   ## Gusmap maps
   out_name <- paste0(SNPCall, "_vcf_gusmap_map.txt")
   map_gus <- create_gusmap_report(vcf_file,"gusmap", "vcf", SNPCall, parent1, parent2)
   write_report(map_gus, out_name)

   out_name <- paste0(SNPCall, "_bam_gusmap_map.txt")
   map_gus <- create_gusmap_report(new.vcf, "gusmap", "bam", SNPCall, parent1, parent2)
   write_report(map_gus, out_name)
   
   RSCRIPT
 >>>

 runtime{
   docker:"taniguti/onemap"
 }

 output{
   File depths_df_vcf_rds = "~{methodName}_df_vcf_depths.rds"
   File depths_updog_vcf_rds = "~{methodName}_updog_vcf_depths.rds"
   File depths_supermassa_vcf_rds = "~{methodName}_supermassa_vcf_depths.rds"
   File depths_polyrad_vcf_rds = "~{methodName}_polyrad_vcf_depths.rds"
   File depths_updog_bam_rds = "~{methodName}_updog_bam_depths.rds"
   File depths_supermassa_bam_rds = "~{methodName}_supermassa_bam_depths.rds"
   File depths_polyrad_bam_rds = "~{methodName}_polyrad_bam_depths.rds"
   File times_vcf_df = "~{methodName}_vcf_df_times.txt"
   File times_vcf_GQ = "~{methodName}_vcf_GQ_times.txt"
   File times_vcf_updog = "~{methodName}_vcf_updog_times.txt"
   File times_vcf_supermassa = "~{methodName}_vcf_supermassa_times.txt"
   File times_vcf_polyrad = "~{methodName}_vcf_polyrad_times.txt"
   File times_vcf_gusmap = "~{methodName}_vcf_gusmap_times.txt"
   File times_bam_gusmap = "~{methodName}_bam_gusmap_times.txt"
   File times_bam_updog = "~{methodName}_bam_updog_times.txt"
   File times_bam_supermassa = "~{methodName}_bam_supermassa_times.txt"
   File times_bam_polyrad = "~{methodName}_bam_polyrad_times.txt"
   File filters_vcf_dfAndGQ = "~{methodName}_filters_vcf_dfAndGQ.txt"
   File filters_vcf_polyrad = "~{methodName}_filters_vcf_polyrad.txt"
   File filters_vcf_supermassa = "~{methodName}_filters_vcf_supermassa.txt"
   File filters_vcf_updog = "~{methodName}_filters_vcf_updog.txt"
   File filters_bam_polyrad = "~{methodName}_filters_bam_polyrad.txt"
   File filters_bam_supermassa = "~{methodName}_filters_bam_supermassa.txt"
   File filters_bam_updog = "~{methodName}_filters_bam_updog.txt"
   File RData_vcf_df = "~{methodName}_vcf_df.RData"
   File RData_vcf_GQ = "~{methodName}_vcf_GQ.RData"
   File RData_vcf_polyrad = "~{methodName}_vcf_polyrad.RData"
   File RData_vcf_supermassa = "~{methodName}_vcf_supermassa.RData"
   File RData_vcf_updog = "~{methodName}_vcf_updog.RData"
   File RData_bam_polyrad = "~{methodName}_bam_polyrad.RData"
   File RData_bam_supermassa = "~{methodName}_bam_supermassa.RData"
   File RData_bam_updog = "~{methodName}_bam_updog.RData"
   File RData_gusmap = "~{methodName}_vcf_gusmap.RData"
   File RData_bam_gusmap = "~{methodName}_bam_gusmap.RData"
   File map_vcf_df = "~{methodName}_vcf_df_map.txt"
   File map_vcf_GQ = "~{methodName}_vcf_GQ_map.txt"
   File map_vcf_polyrad = "~{methodName}_vcf_polyrad_map.txt"
   File map_vcf_supermassa = "~{methodName}_vcf_supermassa_map.txt"
   File map_vcf_updog = "~{methodName}_vcf_updog_map.txt"
   File map_bam_polyrad = "~{methodName}_bam_polyrad_map.txt"              
   File map_bam_supermassa = "~{methodName}_bam_supermassa_map.txt"
   File map_bam_updog = "~{methodName}_bam_updog_map.txt"
   File map_gusmap = "~{methodName}_vcf_gusmap_map.txt"
   File map_bam_gusmap = "~{methodName}_bam_gusmap_map.txt"
 }
}

task JointDatas{
    input{
       Array[File] depths_df_vcf_rds 
       Array[File] depths_updog_vcf_rds 
       Array[File] depths_supermassa_vcf_rds 
       Array[File] depths_polyrad_vcf_rds 
       Array[File] depths_updog_bam_rds 
       Array[File] depths_supermassa_bam_rds 
       Array[File] depths_polyrad_bam_rds 
       Array[File] times_vcf_df 
       Array[File] times_vcf_GQ 
       Array[File] times_vcf_updog 
       Array[File] times_vcf_supermassa 
       Array[File] times_vcf_polyrad 
       Array[File] times_vcf_gusmap 
       Array[File] times_bam_gusmap 
       Array[File] times_bam_updog 
       Array[File] times_bam_supermassa 
       Array[File] times_bam_polyrad 
       Array[File] filters_vcf_dfAndGQ 
       Array[File] filters_vcf_polyrad 
       Array[File] filters_vcf_supermassa 
       Array[File] filters_vcf_updog 
       Array[File] filters_bam_polyrad 
       Array[File] filters_bam_supermassa 
       Array[File] filters_bam_updog 
       Array[File] RData_vcf_df 
       Array[File] RData_vcf_GQ 
       Array[File] RData_vcf_polyrad 
       Array[File] RData_vcf_supermassa 
       Array[File] RData_vcf_updog 
       Array[File] RData_bam_polyrad 
       Array[File] RData_bam_supermassa 
       Array[File] RData_bam_updog 
       Array[File] RData_gusmap 
       Array[File] RData_bam_gusmap 
       Array[File] map_vcf_df 
       Array[File] map_vcf_GQ 
       Array[File] map_vcf_polyrad 
       Array[File] map_vcf_supermassa 
       Array[File] map_vcf_updog 
       Array[File] map_bam_polyrad 
       Array[File] map_bam_supermassa 
       Array[File] map_bam_updog 
       Array[File] map_gusmap 
       Array[File] map_bam_gusmap 
    }

    command <<<

       R --vanilla --no-save <<RSCRIPT
                        
       system("cp ~{sep = " " depths_df_vcf_rds} .")
       system("cp ~{sep = " " depths_updog_vcf_rds } .")
       system("cp ~{sep = " " depths_supermassa_vcf_rds } .")
       system("cp ~{sep = " " depths_polyrad_vcf_rds } .")
       system("cp ~{sep = " " depths_updog_bam_rds } .")
       system("cp ~{sep = " " depths_supermassa_bam_rds } .")
       system("cp ~{sep = " " depths_polyrad_bam_rds } .")
       system("cp ~{sep = " " times_vcf_df } .")
       system("cp ~{sep = " " times_vcf_GQ } .")
       system("cp ~{sep = " " times_vcf_updog } .")
       system("cp ~{sep = " " times_vcf_supermassa } .")
       system("cp ~{sep = " " times_vcf_polyrad } .")
       system("cp ~{sep = " " times_vcf_gusmap } .")
       system("cp ~{sep = " " times_bam_updog } .")
       system("cp ~{sep = " " times_bam_supermassa } .")
       system("cp ~{sep = " " times_bam_polyrad } .")
       system("cp ~{sep = " " times_bam_gusmap } .")
       system("cp ~{sep = " " filters_vcf_dfAndGQ } .")
       system("cp ~{sep = " " filters_vcf_polyrad } .")
       system("cp ~{sep = " " filters_vcf_supermassa } .")
       system("cp ~{sep = " " filters_vcf_updog } .")
       system("cp ~{sep = " " filters_bam_polyrad } .")
       system("cp ~{sep = " " filters_bam_supermassa } .")
       system("cp ~{sep = " " filters_bam_updog } .")
       system("cp ~{sep = " " RData_vcf_df } .")
       system("cp ~{sep = " " RData_vcf_GQ } .")
       system("cp ~{sep = " " RData_vcf_polyrad } .")
       system("cp ~{sep = " " RData_vcf_supermassa } .")
       system("cp ~{sep = " " RData_vcf_updog } .")
       system("cp ~{sep = " " RData_bam_polyrad } .")
       system("cp ~{sep = " " RData_bam_supermassa } .")
       system("cp ~{sep = " " RData_bam_updog } .")
       system("cp ~{sep = " " RData_gusmap } .")
       system("cp ~{sep = " " RData_bam_gusmap } .")
       system("cp ~{sep = " " map_vcf_df } .")
       system("cp ~{sep = " " map_vcf_GQ } .")
       system("cp ~{sep = " " map_vcf_polyrad } .")
       system("cp ~{sep = " " map_vcf_supermassa } .")
       system("cp ~{sep = " " map_vcf_updog } .")
       system("cp ~{sep = " " map_bam_polyrad } .")
       system("cp ~{sep = " " map_bam_supermassa } .")
       system("cp ~{sep = " " map_bam_updog } .")
       system("cp ~{sep = " " map_gusmap } .")
       system("cp ~{sep = " " map_bam_gusmap } .")

       ## Map, times and RDatas
                        
       SNPCall <- c("gatk", "freebayes")
       GenoCall <- c("supermassa", "updog", "df", "polyrad", "GQ", "gusmap")
       CountsFrom <- c("bam", "vcf")
                        
       data_times <- data_map <- temp <- temp.1 <- data.frame()
       z <- 1

       names_list <- vector()
       data_RDatas <- list()
         for(i in SNPCall){
            for(j in GenoCall){
               for(w in CountsFrom){
                 if((j == "df" | j =="GQ") && w=="bam"){
                 } else {
                  temp <- read.table(paste0(i,"_",w,"_",j,"_map.txt"), header=T)
                  temp.1 <- read.table(paste0(i,"_",w,"_",j,"_times.txt"), header=T)
                  load(paste0(i,"_",w,"_",j,".RData"))
                 }
                 data_map <- rbind(data_map, temp)
                 data_times <- rbind(data_times, temp.1)
                 if(j != "gusmap"){
                    data_RDatas[[z]] <- map_out
                 } else {
                    data_RDatas[[z]] <- mydata
                 }
                 names_list <- c(names_list,paste0(i,"_",w,"_",j))
                 z <- z + 1
                }
             }
          }
                         
        saveRDS(data_map, "data_map.rds")
        saveRDS(data_times, "data_times.rds")
        
        names(data_RDatas) <- names_list
        save(data_RDatas, file = "data_RDatas.RData")
        
        # Filters
        
        GenoCall <- c("supermassa", "updog", "dfAndGQ", "polyrad")
        
        data_filters <- temp <- data.frame()
        for(i in SNPCall){
            for(j in GenoCall){
                for(w in CountsFrom){
                    if((j == "dfAndGQ") & w=="bam"){
                    } else {
                     temp <- read.table(paste0(i,"_filters_",w,"_",j,".txt"), header=T)
                    }
                  data_filters <- rbind(data_filters, temp)
                }
            }
        }
        
        saveRDS(data_filters, "data_filters.rds")
        
       # Depths
        
        GenoCall <- c("supermassa", "updog", "df", "polyrad")
        
        data_depths <- data.frame()
        for(i in SNPCall){
            for(j in GenoCall){
                for(w in CountsFrom){
                   if(j == "df" & w == "bam"){
                   } else {
                     temp <- readRDS(paste0(i,"_", j, "_", w, "_depths.rds"))
                     temp <- data.frame(SNPCall=i, GenoCall=i, CountsFrom=w, temp)
                     data_depths <- rbind(data_depths, temp)
                   }
                }
            }
        }
                     
       saveRDS(data_depths, "data_depths.rds")
                         
       RSCRIPT
     >>>

    runtime{
      docker:"taniguti/onemap"
    }

    output{
      File data_depths = "data_depths.rds"
      File data_filters = "data_filters.rds"
      File data_map = "data_map.rds"
      File data_times = "data_times.rds"
      File data_RDatas = "data_RDatas.RData"
    }
}