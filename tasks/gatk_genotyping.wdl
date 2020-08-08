version 1.0

import "../structs/alignment_struct.wdl"
import "../structs/reads_simuS.wdl"
import "../structs/snpcalling_empS.wdl"
import "./utils.wdl" as utils
import "./utilsR.wdl" as utilsR
import "split_filt_vcf.wdl" as norm_filt

workflow GatkGenotyping {
  input {
    Array[Alignment] alignments
    ReferenceFasta references
    String program
    SplitVCF splitvcf
    Array[String] sampleNames
  }

  scatter (alignment in alignments) {
    call HaplotypeCallerERC {
      input:
        ref        = references.ref_fasta,
        geno_fai   = references.ref_fasta_index,
        sampleName = alignment.sample,
        bam_rg     = alignment.bam,
        bam_rg_idx = alignment.bai,
        geno_dict  = references.ref_dict
    }
  }

  call CreateGatkDatabase{
    input:
      path_gatkDatabase = "my_database",
      GVCFs             = HaplotypeCallerERC.GVCF,
      GVCFs_idx         = HaplotypeCallerERC.GVCF_idx,
      ref               = references.ref_fasta
  }

  call GenotypeGVCFs {
    input:
      workspace_tar = CreateGatkDatabase.workspace_tar,
      fasta=references.ref_fasta,
      fasta_fai=references.ref_fasta_index,
      fasta_dict=references.ref_dict
  }
  
  call norm_filt.SplitFiltVCF{
    input:
      vcf_in=GenotypeGVCFs.vcf,
      program=program,
      chromosome = splitvcf.chromosome,
      reference = references.ref_fasta,
      reference_idx = references.ref_fasta_index,
      parent1 = splitvcf.parent1,
      parent2 = splitvcf.parent2
  }

  scatter (alignment in alignments) {
    call utils.BamCounts {
      input:
        sample=alignment.sample,
        program=program,
        bam=alignment.bam,
        bai=alignment.bai,
        ref=references.ref_fasta,
        ref_fai=references.ref_fasta_index,
        ref_dict=references.ref_dict,
        vcf=SplitFiltVCF.vcf_bi_chr_norm,
        tbi=SplitFiltVCF.vcf_bi_chr_norm_tbi
    }
  }
  
  call utils.BamCounts4Onemap {
    input:
      sampleName=sampleNames,
      counts=BamCounts.counts,
      method = program
  }

  call utilsR.BamDepths2Vcf{
    input:
      vcf_file = SplitFiltVCF.vcf_bi_chr_norm,
      ref_bam = BamCounts4Onemap.ref_bam,
      alt_bam = BamCounts4Onemap.alt_bam,
      example_alleles = BamCounts4Onemap.ref_alt_alleles,
      program = program
  }

  output {
    File vcf = SplitFiltVCF.vcf_bi_chr_norm
    File tbi = SplitFiltVCF.vcf_bi_chr_norm_tbi
    File vcf_bi_tot = SplitFiltVCF.vcf_bi_norm
    File vcf_multi_tot = SplitFiltVCF.vcf_multi_norm
    File vcf_bi_bam_counts = BamDepths2Vcf.bam_vcf
    File alt_bam = BamCounts4Onemap.alt_bam
    File ref_bam = BamCounts4Onemap.ref_bam
  }
}


task HaplotypeCallerERC {
  input {
    File ref
    File geno_fai
    String sampleName
    File bam_rg
    File bam_rg_idx
    File geno_dict
  }

  command <<<
    /gatk/gatk HaplotypeCaller \
      -ERC GVCF \
      -R ~{ref} \
      -I ~{bam_rg} \
      -O ~{sampleName}_rawLikelihoods.g.vcf \
      --max-reads-per-alignment-start 0
      
  >>>

  runtime {
    docker: "taniguti/gatk-picard"
    mem:"--nodes=1"
    cpu:1
    time:"120:00:00"
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
    File ref
  }

  command <<<
     
     grep ">" ~{ref} > interval_list_temp
     sed 's/^.//' interval_list_temp > interval.list 
      
     /gatk/gatk GenomicsDBImport \
        --genomicsdb-workspace-path ~{path_gatkDatabase} \
        -L interval.list \
        -V ~{sep=" -V "  GVCFs}

     tar -cf ~{path_gatkDatabase}.tar ~{path_gatkDatabase}
     
  >>>

  runtime {
      docker: "taniguti/gatk-picard"
      mem:"--nodes=1"
      cpu:1
      time:"120:00:00"
  }

  output {
      File workspace_tar = "${path_gatkDatabase}.tar"
  }
}

# Variant calling on gVCF
task GenotypeGVCFs {

  input {
    File workspace_tar
    File fasta
    File fasta_fai
    File fasta_dict
  }

  command <<<
    tar -xf ~{workspace_tar}
    WORKSPACE=$( basename ~{workspace_tar} .tar)
    
    /gatk/gatk GenotypeGVCFs \
        -R ~{fasta} \
        -O gatk.vcf.gz \
        -G StandardAnnotation \
        -V gendb://$WORKSPACE
    
  >>>

  runtime {
    docker: "taniguti/gatk-picard"
    mem:"--nodes=1"
    cpu:1
    time:"120:00:00"
  }

  output {
    File vcf = "gatk.vcf.gz"
    File tbi = "gatk.vcf.gz.tbi"
  }
}
