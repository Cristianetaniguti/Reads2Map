version 1.0

import "../structs/dna_seq_structs.wdl"
import "../tasks/tassel.wdl"
import "../tasks/BWA.wdl"
import "../tasks/utils.wdl" as utils

import "../subworkflows/norm_filt_vcf.wdl" as norm_filt

workflow TasselGenotyping {
    input {
        File families_info
        ReferenceFasta references
        Int max_cores
        Int max_ram
        String enzyme = "ApeKI"
    }

    call tassel.transposeSamples {
        input:
            families_info = families_info
    }

    Array[Array[String]] sample_file = read_tsv(transposeSamples.transpose_samples)

    call tassel.BarcodeFaker {
        input:
            fastq = sample_file[0],
            FullSampleName = sample_file[2]
    }

    call tassel.TasselBeforeAlign {
        input:
            fastq = BarcodeFaker.barcode_fastq,
            enzyme = enzyme,
            key_file = BarcodeFaker.key_file,
            max_ram = max_ram
    }

    call BWA.RunBwaAlignment {
        input:
            sampleName = ["tasselSample"],
            reads1 = [TasselBeforeAlign.fastq_align],
            pair_end = false,
            libraries = ["tasselSample"],
            references = references,
            max_cores = max_cores,
            rm_dupli = false
    }

    call tassel.TasselAfterAlign {
        input:
            tassel_database = TasselBeforeAlign.tassel_database,
            bam = RunBwaAlignment.bam[0],
            max_ram = max_ram,
            enzyme = enzyme,
            key_file = BarcodeFaker.key_file,
            fastq = BarcodeFaker.barcode_fastq
    }

    call norm_filt.Normalization {
        input:
            vcf_in= TasselAfterAlign.tassel_vcf,
            reference = references.ref_fasta,
            reference_idx = references.ref_fasta_index,
            reference_dict = references.ref_dict,
            program = "tassel",
            counts_source = "vcf",
            ploidy = "2"
    }
    
    output {
        Array[File] vcfs = [Normalization.vcf_norm]
        Array[String] software_sele = ["tassel"]
        Array[String] source_sele = ["vcf"]
    }
}