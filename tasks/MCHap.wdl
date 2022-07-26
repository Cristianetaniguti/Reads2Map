version 1.0

workflow MCHap{
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

  call BamToBed {
      input:
        merged_bams = merged_bams
  }

  call SepareChunksBed {
    input:
     bed_file = BamToBed.merged_bed,
     n_nodes = n_nodes
  }

  # If running outside of Reads2Map workflow
  #Array[File] bams = read_lines(bam_list)
  #Array[File] bais = read_lines(bais_list)

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

    call OneMCHap_recall {
        input:
            bams = map_bams["bam"],
            bais = map_bams["bai"],
            vcf_file = OneMCHap.assemble_vcf,
            ploidy = ploidy,
            max_cores = max_cores  
    }
  }

  call mergeVCFs {
      input:
        haplo_vcf = OneMCHap_recall.haplo_vcf
  }

  call FilterMulti {
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

task BamToBed {
    input {
        File? merged_bams
    }

    command <<<
        bamToBed -i ~{merged_bams} > file.bed
        bedtools merge -i file.bed > merged.bed
    >>>

    runtime {
        job_name: "BamToBed"
        docker: "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
        node:"--nodes=1"
        mem:"--mem=10G"
        tasks:"--ntasks=1"
        time:"05:00:00"
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Uses [bamToBed](https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html) and [bedtools](https://bedtools.readthedocs.io/en/latest/) to create BED file and merge overlapping intervals."
    }    

    output {
        File merged_bed = "merged.bed"
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

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Split BED file rows in batches according to defined number of nodes."
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
            --mcmc-steps 2000 \
            --haplotype-posterior-threshold 0.9 \
            --cores ~{max_cores} | bgzip > assemble.vcf.gz
            
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

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Runs [MCHap](https://github.com/PlantandFoodResearch/MCHap) assemble step in a batch of the BED file rows."
    }  

    output {
        File assemble_vcf = "assemble.vcf.gz"
    }
}

task OneMCHap_recall {
    input{
        Array[File] bams
        Array[File] bais
        File vcf_file
        Int ploidy
        Int max_cores
    }

    command <<<

        export TMPDIR=/tmp

        ln -s ~{sep=" " bams} .
        ln -s ~{sep=" " bais} .

        tabix -p vcf ~{vcf_file}
            
        mchap call \
            --haplotypes ~{vcf_file} \
            --bam *.bam \
            --ploidy ~{ploidy} \
            --inbreeding 0.01 \
            --base-error-rate 0.0025 \
            --ignore-base-phred-scores \
            --mcmc-burn 1000 \
            --mcmc-steps 2000 \
            --cores ~{max_cores} \
            | bgzip > haplotypes.vcf.gz

    >>>

    runtime {
        docker: "cristaniguti/mchap:0.0.1"
        # memory: "4 GB"
        # cpu: 1
        # preemptible: 3
        # disks: "local-disk " + disk_size + " HDD"
        job_name: "MCHap_recall"
        node:"--nodes=1"
        mem:"--mem=20G" 
        tasks:"--ntasks-per-node=16"
        time:"24:00:00"
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Runs [MCHap](https://github.com/PlantandFoodResearch/MCHap) call step in a batch of the BED file rows."
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

        ln -s ~{sep=" " haplo_vcf} .

        for file in $(echo haplotypes*); do
            filename=$(basename -- "$file")
            name="${filename%.*}"
            echo $name
            bcftools sort $file --output-file $name.sorted.vcf
            bgzip $name.sorted.vcf
            tabix -p vcf $name.sorted.vcf.gz
        done

        bcftools concat *sorted.vcf.gz --output merged.vcf.gz
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

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Uses [bcftools](https://samtools.github.io/bcftools/bcftools.html) to sort and merge VCF files."
    }

    output {
        File merged_vcf = "merged.sorted.vcf.gz"
    }
}

task FilterMulti {
    input {
        File multi_vcf
        String? P1
        String? P2
        Int ploidy
    }

    command <<<
        R --vanilla --no-save <<RSCRIPT

            library(Reads2MapTools)
            filter_multi_vcf("~{multi_vcf}", "~{P1}", "~{P2}", 
                             ploidy = ~{ploidy},
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

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Filters VCF file markers according to segregation expected in a outcrossing F1 population. Adapts the alleles codification. See [Reads2MapTools](https://github.com/Cristianetaniguti/Reads2MapTools) for more information."
    }

    output {
        File multi_vcf_filt = "multi_vcf_filt.vcf.gz"
    }
}