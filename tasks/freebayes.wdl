version 1.0

task RunFreebayes {

  input {
    File reference
    File reference_idx
    Array[File] bam
    Array[File] bai
    Int max_cores
  }

  Int disk_size = ceil(size(reference, "GiB") + size(bam, "GiB") +  50)
  Int memory_size = ceil(size(bam, "MiB") * 1.25)

  command <<<
   # needed for some singularity versions
   ln -sf ~{sep=" " bam} .
   ln -sf ~{sep=" " bai} .

   freebayes-parallel <(fasta_generate_regions.py ~{reference_idx} 100000) ~{max_cores} \
   --genotype-qualities -f ~{reference} *.bam > "freebayes.vcf"

  >>>

  runtime {
    docker: "cristaniguti/freebayes:0.0.1"
    cpu: max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "RunFreebayes"
    mem:"~{memory_size}M"
    time:"48:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Split genomic regions and runs [freebayes](https://github.com/freebayes/freebayes) parallelized."
  }

  output {
    File vcf = "freebayes.vcf"
  }
}
