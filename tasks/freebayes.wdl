version 1.0

task RunFreebayes {

  input {
    File reference
    File reference_idx
    File bam
    File bai
    Int max_cores
    Int ploidy
  }

  Int disk_size = ceil(size(reference, "GiB") + size(bam, "GiB") +  50)
  Int memory_size = ceil(size(bam, "MiB") * 25 * max_cores + 10000)

  command <<<

   ln -s ~{bam} .
   ln -s ~{bai} .

   freebayes-parallel <(fasta_generate_regions.py ~{reference_idx} 100000) ~{max_cores} \
   --genotype-qualities --ploidy ~{ploidy} -f ~{reference} *bam > "freebayes.vcf"

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
