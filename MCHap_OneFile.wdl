version 1.0

import "tasks/MCHap.wdl" as mchap

workflow MCHap{
  input {
    File reference
    File reference_idx
    File vcf_file
    Int n_nodes
    Int max_cores
    File bam_list 
    File bais_list 
    File merged_bed 
    Int ploidy
    String P1
    String P2
  }

  call mchap.SepareChunksBed {
    input:
     bed_file = merged_bed,
     n_nodes = n_nodes
  }

  # If running outside of Reads2Map workflow
  Array[File] bams = read_lines(bam_list)
  Array[File] bais = read_lines(bais_list)

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

  call mchap.mergeVCFs {
      input:
        haplo_vcf = OneMCHap_recall.haplo_vcf
  }

  call mchap.FilterMulti {
      input:
        multi_vcf = mergeVCFs.merged_vcf,
        ploidy = ploidy,
        P1 = P1,
        P2 = P2
  }

  output {
    File haplo_vcf_merged = FilterMulti.multi_vcf_filt
  }
}


