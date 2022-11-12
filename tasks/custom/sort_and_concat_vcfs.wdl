version 1.0

task mergeVCFs {
    input {
        Array[File] haplo_vcf
    }

    Int disk_size = ceil(size(haplo_vcf, "GiB") * 1.5)
    Int memory_size = 5000

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
        cpu: 1
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "mergeVCFs"
        mem:"~{memory_size}M"
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
