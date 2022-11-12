version 1.0


# It will always produce P1, P2, F1 and then F2_00X, where
# X will increase from 1 to samples
task GenerateSampleNames {  # TODO: probably a name like 'ReadSamplesNamesInVcf' is better

  input {
    File simulated_vcf
  }

  Int disk_size = ceil(size(simulated_vcf, "GiB") * 2) 
  Int memory_size = 1000 

  command <<<
    export PATH=$PATH:/opt/conda/bin

    python <<CODE
    from pysam import VariantFile

    bcf_in = VariantFile("~{simulated_vcf}")

    for i in bcf_in.header.samples:
        print(i)
    CODE

  >>>

  runtime {
    docker: "cristaniguti/miniconda-alpine:0.0.1"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "GenerateSampleNames"
    mem:"~{memory_size}M"
    time:"05:00:00"
  }

  meta {
    author: "Lucas Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Creates the sample names."
  }

  output {
    Array[String] names = read_lines(stdout())
  }
}
