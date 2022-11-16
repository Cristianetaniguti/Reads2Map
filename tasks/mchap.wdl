version 1.0

task OneMCHap {
    input {
        Array[File] bams
        Array[File] bais
        File bed
        File vcf_file
        File reference
        File reference_idx
        Int ploidy
        Int max_cores
    }

    Int disk_size = ceil(size(bams, "GiB") * 1.5 + size(bed, "GiB") * 1.5 + size(vcf_file, "GiB") * 1.5 + size(reference, "GiB"))
    Int memory_size = 3000

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
        cpu: max_cores
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "MCHap"
        mem:"~{memory_size}M"
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
    input {
        Array[File] bams
        Array[File] bais
        File vcf_file
        Int ploidy
        Int max_cores
    }

    Int disk_size = ceil(size(bams, "GiB") * 1.25 + size(vcf_file, "GiB") * 1.5)
    Int memory_size = 3000

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
        cpu: max_cores
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "OneMCHap_recall"
        mem:"~{memory_size}M"
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
