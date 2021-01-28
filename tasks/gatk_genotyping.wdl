version 1.0

import "../structs/snpcalling_empS.wdl"
import "../structs/reference_struct.wdl"
import "./utils.wdl" as utils
import "./utilsR.wdl" as utilsR
import "split_filt_vcf.wdl" as norm_filt

workflow GatkGenotyping {
  input {
    Array[File] bam
    Array[File] bai
    Reference references
    String program
    String parent1
    String parent2
    String chrom
    Array[String] sampleNames
    Int max_cores
  }

  call HaplotypeCallerERC {
     input:
       ref        = references.ref_fasta,
       geno_fai   = references.ref_fasta_index,
       bam_rg     = bam,
       geno_dict  = references.ref_dict
   }

  Map[String, Array[File]] vcfs = {"vcf": HaplotypeCallerERC.GVCF, "idx": HaplotypeCallerERC.GVCF_idx}

  call CreateGatkDatabase{
    input:
      path_gatkDatabase = "my_database",
      GVCFs             = vcfs["vcf"],
      GVCFs_idx         = vcfs["idx"],
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
      reference = references.ref_fasta,
      reference_idx = references.ref_fasta_index,
      parent1 = parent1,
      parent2 = parent2
  }

  Map[String, Array[File]] bams = {"bam": bam, "bai": bai}

  call utils.ReplaceAD {
    input:
      ref_fasta = references.ref_fasta,
      ref_index = references.ref_fasta_index,
      bams = bams["bam"],
      bais = bams["bai"],
      vcf = SplitFiltVCF.vcf_bi,
      tbi = SplitFiltVCF.vcf_bi_tbi,
      program = program
  }

  output {
    File vcf_bi = SplitFiltVCF.vcf_bi
    File tbi_bi = SplitFiltVCF.vcf_bi_tbi
    File vcf_multi = SplitFiltVCF.vcf_multi
    File vcf_bi_bam_counts = ReplaceAD.bam_vcf
  }
}


task HaplotypeCallerERC {
  input {
    File ref
    File geno_fai
    Array[File] bam_rg
    File geno_dict
  }

  command <<<

    for bam in ~{sep= " " bam_rg}; do

      samtools index $bam
      sample=`basename -s .sorted.bam $bam`
      echo $sample

      /gatk/gatk HaplotypeCaller \
        -ERC GVCF \
        -R ~{ref} \
        -I "$bam" \
        -O "rawLikelihoods.g.vcf" \
        --max-reads-per-alignment-start 0

      mv rawLikelihoods.g.vcf $sample.rawLikelihoods.g.vcf
      mv rawLikelihoods.g.vcf.idx $sample.rawLikelihoods.g.vcf.idx

    done

  >>>

  runtime {
    docker: "taniguti/gatk-picard"
    mem:"20GB"
    cpu:1
    time:"14:00:00"
  }

  output {
    Array[File] GVCF = glob("*.rawLikelihoods.g.vcf")
    Array[File] GVCF_idx = glob("*.rawLikelihoods.g.vcf.idx")
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

     ln -sf ~{sep=" " GVCFs} .
     ln -sf ~{sep=" " GVCFs_idx} .
     
     VCFS=$(echo *.g.vcf)
     VCFS=${VCFS// / -V }

     /gatk/gatk GenomicsDBImport \
        --genomicsdb-workspace-path ~{path_gatkDatabase} \
        -L interval.list \
        -V $VCFS

     tar -cf ~{path_gatkDatabase}.tar ~{path_gatkDatabase}

  >>>

  runtime {
      docker: "taniguti/gatk-picard"
      mem:"30GB"
      cpu:1
      time:"15:00:00"
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
    mem:"30GB"
    cpu:1
    time:"10:00:00"
  }

  output {
    File vcf = "gatk.vcf.gz"
    File tbi = "gatk.vcf.gz.tbi"
  }
}
