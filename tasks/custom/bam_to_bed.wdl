version 1.0

task BamToBed {
    input {
        File? merged_bams
    }

    Int disk_size = ceil(size(merged_bams, "GiB") * 1.5)
    Int memory_size = 3000

    command <<<
        bamToBed -i ~{merged_bams} > file.bed
        bedtools merge -i file.bed > merged.bed
    >>>

    runtime {
        docker: "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
        cpu: 1
        # Cloud
        memory:"~{memory_size} MiB"
        disks:"local-disk " + disk_size + " HDD"
        # Slurm
        job_name: "BamToBed"
        mem:"~{memory_size}M"
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
