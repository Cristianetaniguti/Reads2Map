version 1.0

workflow MCHap{
  input {
    File reference
    File reference_idx
    File vcf_file
    File vcf_tbi
    File bed_file
    Int n_nodes
    Int max_cores
    File bam_list
    File bais_list
    Int ploidy
  }

  call SepareChunksBed {
    input:
     bed_file = bed_file,
     n_nodes = n_nodes
  }

  Array[File] bams = read_lines(bam_list)
  Array[File] bais = read_lines(bais_list)

  scatter (bed_chunk in SepareChunksBed.chunks){
    call OneMCHap {
        input:
        bams = bams,
        bais = bais,
        bed = bed_chunk,
        vcf_file = vcf_file,
        vcf_tbi = vcf_tbi,
        reference = reference,
        reference_idx = reference_idx,
        ploidy = ploidy,
        max_cores = max_cores
    }
  }

  call mergeVCFs {
      input:
        haplo_vcf = OneMCHap.haplo_vcf
  }

  output {
    File haplo_vcf_merged = mergeVCFs.merged_vcf
  }
}

task SepareChunksBed {
    input {
        File bed_file
        Int n_nodes
    }

    command <<<
        R --vanilla --no-save <<RSCRIPT
            library(vroom)
            bed <- vroom("~{bed_file}")
            
            chunk_size <- floor(dim(bed)[1]/~{n_nodes})
            
            chunk_temp <- rep(1:~{n_nodes}, each = chunk_size)
            chunk <- c(chunk_temp, rep(~{n_nodes}+1, dim(bed)[1] - length(chunk_temp)))

            chunk_sep <- split.data.frame(bed, chunk)

            for(i in 1:length(chunk_sep)){
                write.table(chunk_sep[[i]], file = paste0("chunk_",i, ".txt"), quote = F, col.names = F, row.names = F, sep="\t")
            }
        RSCRIPT
    >>>

    runtime {
        job_name: "SepareChunksBed"
        docker: "cristaniguti/reads2map:0.0.1"
        node:"--nodes=1"
        mem:"--mem=1G"
        tasks:"--ntasks=1"
        time:"00:05:00"
    }

    output {
        Array[File] chunks = glob("chunk*")
    }
}

task OneMCHap {
    input{
        Array[File] bams
        Array[File] bais
        File bed
        File vcf_file
        File vcf_tbi
        File reference
        File reference_idx
        Int ploidy
        Int max_cores
    }

    command <<<

        export TMPDIR=/tmp

        mv ~{sep=" " bams} .
        mv ~{sep=" " bais} .
        mv ~{vcf_file} .
        mv ~{vcf_tbi} .
        mv ~{reference} .
        mv ~{reference_idx} .

        referenceName=$(basename ~{reference})
        vcfName=$(basename ~{vcf_file})

        mchap assemble \
            --bam *.bam \
            --targets ~{bed} \
            --variants $vcfName \
            --reference $referenceName \
            --ploidy ~{ploidy} \
            --cores ~{max_cores} | bgzip > haplotypes.vcf.gz
            
    >>>

    runtime {
        docker: "cristaniguti/mchap:0.0.1"
        # memory: "4 GB"
        # cpu: 1
        # preemptible: 3
        # disks: "local-disk " + disk_size + " HDD"
        job_name: "MCHap"
        node:"--nodes=1"
        mem:"--mem=20G" 
        tasks:"--ntasks-per-node=16"
        time:"24:00:00"
    }

    output {
        File haplo_vcf = "haplotypes.vcf.gz"
    }
}

task mergeVCFs {
    input {
        Array[File] haplo_vcf
    }

    command <<<

        bcftools concat ~{sep=" " haplo_vcf} --output merged.vcf.gz
        bcftools sort merged.vcf.gz --output merged.sorted.vcf.gz

    >>>

    runtime {
        docker:"lifebitai/bcftools:1.10.2"
        # memory: "2 GB"
        # cpu:1
        # preemptible: 3
        job_name: "mergeVCFs"
        node:"--nodes=1"
        mem:"--mem=15GB"
        tasks:"--ntasks=1"
        time:"01:00:00"
    }

    output {
        File merged_vcf = "merged.sorted.vcf.gz"
    }
}