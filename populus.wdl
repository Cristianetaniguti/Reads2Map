version 1.0

import "./structs/populusS.wdl"

workflow populus{

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
  
    bwa mem -t 5 ~{ref} ~{reads1} | \
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

