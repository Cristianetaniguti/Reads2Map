version 1.0

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

  call SepareChunksBed {
    input:
     bed_file = merged_bed,
     n_nodes = n_nodes
  }

  # If running outside of Reads2Map workflow
  Array[File] bams = read_lines(bam_list)
  Array[File] bais = read_lines(bais_list)

  Map[String, Array[File]] map_bams = {"bam": bams, "bai": bais}

  scatter (bed_chunk in SepareChunksBed.chunks){
    call OneMCHap {
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
  }

  call mergeVCFs {
      input:
        haplo_vcf = OneMCHap.haplo_vcf
  }

  call FilterMulti {
      input:
        multi_vcf = mergeVCFs.merged_vcf,
        P1 = P1,
        P2 = P2
  }

  output {
    File haplo_vcf_merged = FilterMulti.multi_vcf_filt
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
        File reference
        File reference_idx
        Int ploidy
        Int max_cores
    }

    command <<<

        export TMPDIR=/tmp

        ln -s ~{reference} .
        ln -s ~{reference_idx} .
        ln -s ~{sep=" " bams} .
        ln -s ~{sep=" " bais} .

        referenceName=$(basename ~{reference})

        tabix -p vcf ~{vcf_file}

        mchap assemble \
            --bam *.bam \
            --targets ~{bed} \
            --variants ~{vcf_file} \
            --reference $referenceName \
            --ploidy ~{ploidy} \
            --inbreeding 0.01 \
            --base-error-rate 0.0025 \
            --ignore-base-phred-scores \
            --mcmc-burn 1000 \
            --mcmc-steps \
            --haplotype-posterior-threshold 0.9 \
            --cores ~{max_cores} | bgzip > assemble.vcf.gz
            
        mchap call
            --haplotypes assemble.vcf.gz \
            --bam *.bam \
            --ploidy ~{ploidy} \
            --inbreeding 0.01 \
            --base-error-rate 0.0025 \
            --ignore-base-phred-scores \
            --mcmc-burn 1000 \
            --mcmc-steps 2000 \
            | bgzip > haplotypes.vcf.gz

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
        bcftools sort merged.vcf.gz --output-file merged.sorted.vcf.gz

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

task FilterMulti {
    input {
        File multi_vcf
        String P1
        String P2
    }

    command <<<
        R --vanilla --no-save <<RSCRIPT

            library(Reads2MapTools)
            filter_multi_vcf("~{multi_vcf}", "~{P1}", "~{P2}", 
                             vcf.out = "multi_vcf_filt.vcf.gz")

        RSCRIPT
    >>>

    runtime {
        docker:"cristaniguti/reads2map:0.0.1"
        # memory: "2 GB"
        # cpu:1
        # preemptible: 3
        job_name: "FilterMulti"
        node:"--nodes=1"
        mem:"--mem=15GB"
        tasks:"--ntasks=1"
        time:"01:00:00"
    }

    output {
        File multi_vcf_filt = "multi_vcf_filt.vcf.gz"
    }
}