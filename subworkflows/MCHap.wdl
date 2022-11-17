version 1.0

import "../tasks/chunk_lists.wdl"
import "../tasks/utilsR.wdl"
import "../tasks/utils.wdl"
import "../tasks/mchap.wdl"

workflow MCHap {
  input {
    File reference
    File reference_idx
    File vcf_file
    Int n_nodes
    Int max_cores
    Array[File] bams # if file change to bam_list
    Array[File] bais # if file change to bais_list
    File? merged_bams
    Int ploidy
    String? P1
    String? P2
  }

  call utils.BamToBed {
      input:
        merged_bams = merged_bams
  }

  call chunk_lists.SepareChunksBed {
    input:
     bed_file = BamToBed.merged_bed,
     n_nodes = n_nodes
  }

  # If running outside of Reads2Map workflow
  #Array[File] bams = read_lines(bam_list)
  #Array[File] bais = read_lines(bais_list)

  Map[String, Array[File]] map_bams = {"bam": bams, "bai": bais}

  scatter (bed_chunk in SepareChunksBed.chunks){
    call mchap.OneMCHap {
        input:
        bams = map_bams["bam"],
        bais = map_bams["bai"],
        bed = bed_chunk,
        vcf_file = vcf_file,
        reference = reference,
        reference_idx = reference_idx,
        ploidy = ploidy,
        max_cores = max_cores
    }

    call mchap.OneMCHap_recall {
        input:
            bams = map_bams["bam"],
            bais = map_bams["bai"],
            vcf_file = OneMCHap.assemble_vcf,
            ploidy = ploidy,
            max_cores = max_cores
    }
  }

  call utils.mergeVCFs {
      input:
        haplo_vcf = OneMCHap_recall.haplo_vcf
  }

   if(defined(P1)){
      call utilsR.FilterMulti {
          input:
            multi_vcf = mergeVCFs.merged_vcf,
            ploidy = ploidy,
            P1 = P1,
            P2 = P2
      }
  }

  File final_vcf = select_first([FilterMulti.multi_vcf_filt, mergeVCFs.merged_vcf])

  output {
    File haplo_vcf_merged = final_vcf
  }
}
