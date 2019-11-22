version 1.0

import "structs/populusS.wdl"

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
      sampleName       = GenerateSampleNames.names,
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
        gatk_alt_depth = BamCounts4Onemap.gatk_alt_bam
    }
  }

  output{
    File filters_dfAndGQ = avalVCFs.filters_dfAndGQ
    File filters_polyrad = avalVCFs.filters_polyrad
    File filters_supermassa = avalVCFs.filters_supermassa
    File filters_updog = avalVCFs.filters_updog
    File filters_bam_polyrad = avalVCFs.filters_bam_polyrad
    File filters_bam_supermassa = avalVCFs.filters_bam_supermassa
    File filters_bam_updog = avalVCFs.filters_bam_updog
    File map_df = avalVCFs.map_df
    File map_GQ = avalVCFs.map_GQ
    File map_polyrad = avalVCFs.map_polyrad
    File map_supermassa = avalVCFs.map_supermassa
    File map_updog = avalVCFs.map_updog
    File map_bam_polyrad = avalVCFs.map_bam_polyrad
    File map_bam_supermassa = avalVCFs.map_bam_supermassa
    File map_bam_updog = avalVCFs.map_bam_updog
    File map_gusmap = avalVCFs.map_gusmap
    File map_bam_gusmap = avalVCFs.map_bam_gusmap
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
    File sampleName
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(R.utils)
      system("cp ~{sep=" "  freebayes_counts} .")
      system("cp ~{sep=" "  gatk_counts} .")
      names <- c("~{sep=" , "  sampleName}")
      names <- unlist(strsplit(names, split = " , "))

      methods <- c("freebayes", "gatk")

      for(method in methods){

      file.counts <- read.table(paste0(names[1],"_", method,"_counts.tsv"), skip = 1448, header=T, stringsAsFactors = F)

      ref_depth_matrix2 <- alt_depth_matrix2  <- matrix(NA, nrow = dim(file.counts)[1], ncol = length(names))

      for(j in 1:length(names)){
        ## From picard tool

        file.counts <- read.table(paste0(names[j],"_", method,"_counts.tsv"), skip = 1448, header=T, stringsAsFactors = F)

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
  }

  command <<<

    R --vanilla --no-save <<RSCRIPT

    system("cp ~{freebayes_ref_depth} .")
    system("cp ~{freebayes_alt_depth} .")
    system("cp ~{gatk_ref_depth} .")
    system("cp ~{gatk_alt_depth} .")
    source("/opt/scripts/functions_empirical.R")

    vcf_file <- "~{vcf_file}"
    method_name <- "~{methodName}"
    ## READING FINAL VCF FROM PIPELINE
    vcf <- read.vcfR(vcf_file)
    df <- onemap_read_vcfR(vcfR.object=vcf,
                          cross= cross,
                          parent1="PT_F",
                          parent2="PT_M")


    ## FILTERS REPORT
    out_name <- paste0(method_name, "_filters_dfAndGQ.txt")
    filters_tab <- create_filters_report(df)
    write_report(filters_tab[[1]], out_name)

    max.cores <- 10 #detectCores() # Change here to the maximun cores available for the process
    n.mk <- length(filters_tab[[2]][[1]])
    min.group.size <- 60

    while(n.mk/max.cores < min.group.size){
      max.cores <- max.cores - 1
    }

    out_name <- paste0(method_name, "_map_df.RData")
    map_out <- parmap(input.seq = filters_tab[[2]], cores = max.cores, overlap = 10)
    save(map_out, out_name)

    # MAPS REPORT - GQ
    aval.gq <- extract_depth(vcfR.object=vcf,
                            onemap.object=df,
                            vcf.par="GQ",
                            parent1="PT_F",
                            parent2="PT_M",
                            f1 = f1,
                            recovering=FALSE)

    aval.gq <- create_probs(df, genotypes_errors=aval.gq)
    filters_tab <- create_filters_report(aval.gq)

    out_name <- paste0(method_name, "_map_GQ.RData")
    map_out <- parmap(input.seq = filters_tab[[2]], cores = max.cores, overlap = 10)
    save(map_out, out_name)

    # OTHER TOOLS
    ## With depths from vcf

    updog.aval <- updog_error(
      vcfR.object=vcf,
      onemap.object=df,
      vcf.par="AD",
      parent1="PT_F",
      parent2="PT_M",
      f1 = f1,
      recovering=TRUE,
      mean_phred=20,
      cores=max.cores,
      depths=NULL)

    supermassa.aval <- supermassa4onemap::supermassa_error(
      vcfR.object=vcf,
      onemap.object = df,
      vcf.par = "AD",
      parent1 = "PT_F",
      parent2 = "PT_M",
      f1 = f1,
      recovering = TRUE,
      mean_phred = 20,
      cores = max.cores,
      depths = NULL)

    polyrad.aval <- polyRAD_error(
      vcf=vcf_file,
      onemap.obj=df,
      parent1="PT_F",
      parent2="PT_M",
      f1 = f1,
      crosstype=cross)

    metodologies <- list(updog = updog.aval, supermassa = supermassa.aval, polyrad = polyrad.aval)
    for (metodology in names(metodologies)){
      error_aval <- metodologies[[metodology]]
      ## Filters
      out_name <- paste0(method_name, "_filters_", metodology, ".txt")
      filters_tab <- create_filters_report(error_aval)
      write_report(filters_tab[[1]], out_name)

      ## Maps
      out_name <- paste0(method_name, "_map_", metodology, ".RData")
      map_out <- parmap(input.seq = filters_tab[[2]], cores = max.cores, overlap = 10)
      save(map_out, out_name)
    }

    ## Depths from bam
    depths.alt <- read.table(paste0(method_name, "_alt_depth_bam.txt"), header = T)
    depths.ref <- read.table(paste0(method_name, "_ref_depth_bam.txt"), header = T)

    depths <- list("ref"=depths.ref, "alt"=depths.alt)

    updog.aval.bam <- updog_error(
      vcfR.object=vcf,
      onemap.object=df,
      vcf.par="AD",
      parent1="PT_F",
      parent2="PT_M",
      f1 = f1,
      recovering=TRUE,
      mean_phred=20,
      cores=max.cores,
      depths=depths)

    supermassa.aval.bam <- supermassa_error(
      vcfR.object=vcf,
      onemap.object = df,
      vcf.par = "AD",
      parent1 = "PT_F",
      parent2 = "PT_M",
      f1 = f1,
      recovering = TRUE,
      mean_phred = 20,
      cores = max.cores,
      depths = depths)

    new.vcf <- make_vcf(vcf_file, depths, method_name)

    polyrad.aval.bam <- polyRAD_error(
      vcf=new.vcf,
      onemap.obj=df,
      parent1="PT_F",
      parent2="PT_M",
      f1 = f1,
      crosstype=cross)

    metodologies <- list(updog = updog.aval.bam, supermassa= supermassa.aval.bam, polyrad=polyrad.aval.bam)
    for (metodology in names(metodologies)){
      error_aval <- metodologies[[metodology]]
      ## Filters
      out_name <- paste0(method_name, "_filters_bam_", metodology, ".txt")
      filters_tab <- create_filters_report(error_aval)
      write_report(filters_tab, out_name)

      ## Maps
      out_name <- paste0(method_name, "_map_bam_", metodology, ".RData")
      map_out <- parmap(input.seq = filters_tab[[2]], cores = max.cores, overlap = 10)
      save(map_out, out_name)
    }

    ## Gusmap maps
    out_name <- paste0(method_name, "_map_gusmap.txt")
    map_gus <- create_gusmap_report(vcf_file)
    write_report(map_gus, out_name)

    out_name <- paste0(method_name, "_map_bam_gusmap.txt")
    map_gus <- create_gusmap_report(new.vcf)
    write_report(map_gus, out_name)

    RSCRIPT

  >>>

  runtime{
    docker:"taniguti/onemap"
  }

  output{
    File filters_dfAndGQ = "~{methodName}_filters_dfAndGQ.txt"
    File filters_polyrad = "~{methodName}_filters_polyrad.txt"
    File filters_supermassa = "~{methodName}_filters_supermassa.txt"
    File filters_updog = "~{methodName}_filters_updog.txt"
    File filters_bam_polyrad = "~{methodName}_filters_bam_polyrad.txt"
    File filters_bam_supermassa = "~{methodName}_filters_bam_supermassa.txt"
    File filters_bam_updog = "~{methodName}_filters_bam_updog.txt"
    File map_df = "~{methodName}_map_df.RData"
    File map_GQ = "~{methodName}_map_GQ.RData"
    File map_polyrad = "~{methodName}_map_polyrad.RData"
    File map_supermassa = "~{methodName}_map_supermassa.RData"
    File map_updog = "~{methodName}_map_updog.RData"
    File map_bam_polyrad = "~{methodName}_map_bam_polyrad.RData"
    File map_bam_supermassa = "~{methodName}_map_bam_supermassa.txt"
    File map_bam_updog = "~{methodName}_map_bam_updog.txt"
    File map_gusmap = "~{methodName}_map_gusmap.txt"
    File map_bam_gusmap = "~{methodName}_map_bam_gusmap.txt"
  }
}
