version 1.0

# Insert into a fasta sequence the variants present in a VCF file
# TODO: Probably this task is not beeing used. Consider removing
task RunVcf2diploid {
  input {
    String sampleName
    File ref_genome
    File simu_vcf
    String chromosome
  }

  Int disk_size = ceil(size(ref_genome, "GiB") * 2 + size(simu_vcf, "GiB")) 
  Int memory_size = 5000  

  command <<<
    java -jar /usr/jars/vcf2diploid.jar -id ~{sampleName} -chr ~{ref_genome} -vcf ~{simu_vcf}

  >>>

  runtime {
    docker: "cristaniguti/java-in-the-cloud:0.0.1"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "RunVcf2diploid"
    mem:"~{memory_size}M"
    time:"05:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Uses [vcf2diploid](https://github.com/abyzovlab/vcf2diploid) to include the VCF genotypes into the FASTA file."
  }

  output {
    File maternal_genomes = "~{chromosome}_${sampleName}_maternal.fa"
    File paternal_genomes = "~{chromosome}_${sampleName}_paternal.fa"
  }

}