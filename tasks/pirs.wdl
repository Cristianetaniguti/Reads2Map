version 1.0

# Creates homologous genome with some variation
# specified with -s and -d
task GenerateAlternativeGenome {
  input {
    Int seed
    File ref_genome
  }

  Int disk_size = ceil(size(ref_genome, "GiB") * 2) 

  command <<<
    /pirs/src/pirs/pirs diploid ~{ref_genome} -s 0.0133 -d 0.0022 -v 0 -o alt --random-seed ~{seed}
  >>>

  runtime {
    docker: "cristaniguti/pirs-ddrad-cutadapt:0.0.1"
    cpu:1
    # Cloud
    memory:"3000 MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "GenerateAlternativeGenome"
    mem:"3000M"
    time:"01:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Uses [pirs](https://github.com/galaxy001/pirs) to create homologous genome with variations when reference VCF is not provided."
  }

  output {
    File alt_fasta = "alt.snp.indel.fa"
    File indels = "alt.indel.lst"
    File snps = "alt.snp.lst"
  }
}
