version 1.0

import "../structs/dna_seq_structs.wdl"

import "../tasks/custom/alignment.wdl"
import "../tasks/custom/chunk_lists.wdl"
import "../tasks/gatk.wdl"


workflow GATK_poly {

  input {
    File samples_info
    ReferenceFasta references
    Int max_cores
    Int chunk_size
    Int ploidy
  }

    call chunk_lists.SepareChunksFastqString as SepareChunks {
        input:
            families_info=samples_info,
            chunk_size = chunk_size
    }

    scatter (chunk in SepareChunks.chunks) {

        Array[Array[String]] sample_file = read_tsv(chunk)

        call alignment.RunBwaAlignmentForGatkPolyploid as RunBwaAlignment {
            input:
                sampleName  = sample_file[1],
                reads       = sample_file[0],
                libraries   = sample_file[2],
                references  = references,
                max_cores   = max_cores,
        }
    }

call chunk_lists.CreateChunksBam as CreateChunks {
    input:
      bams=flatten(RunBwaAlignment.bam),
      bams_index=flatten(RunBwaAlignment.bai),
      chunk_size=chunk_size,
      reference_fasta = references.ref_fasta
  }

  scatter (chunk in zip(CreateChunks.bams_chunks, CreateChunks.bais_chunks)) {

    call gatk.HaplotypeCaller as HaplotypeCaller {
      input:
        bams = read_lines(chunk.left),
        bams_index = read_lines(chunk.right),
        reference_fasta = references.ref_fasta,
        reference_fai = references.ref_fasta_index,
        reference_dict = references.ref_dict,
        ploidy = ploidy,
        chunk_size = chunk_size
    }
  }

  Array[String] calling_intervals = read_lines(CreateChunks.interval_list)

  scatter (interval in calling_intervals) {
    call gatk.ImportGVCFs as ImportGVCFs {
      input:
        vcfs=flatten(HaplotypeCaller.vcfs),
        vcfs_index=flatten(HaplotypeCaller.vcfs_index),
        reference_fasta=references.ref_fasta,
        reference_fai=references.ref_fasta_index,
        reference_dict=references.ref_dict,
        interval = interval
    }

    call gatk.GenotypeGVCFs as GenotypeGVCFs {
      input:
        workspace_tar = ImportGVCFs.output_workspace,
        interval = interval,
        reference_fasta = references.ref_fasta,
        reference_fai = references.ref_fasta_index,
        reference_dict = references.ref_dict
    }
  }

  call gatk.MergeVCFs as MergeVCFs {
    input:
      input_vcfs = GenotypeGVCFs.vcf,
      input_vcf_indices = GenotypeGVCFs.vcf_tbi
  }

  output {
    File gatk_vcf = MergeVCFs.output_vcf
  }
}
