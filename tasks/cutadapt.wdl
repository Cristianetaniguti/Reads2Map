version 1.0

task RemoveAdapt {
  input{
    File sequence
    String adapter
    String sequence_name
  }

  command <<<
    cutadapt -a ~{adapter} -o ~{sequence_name}_trim.fastq.gz ~{sequence} --minimum-length 64
  >>>

  runtime {
    docker:"kfdrc/cutadapt"
    job_name: "RemoveAdapt"
    node:"--nodes=1"
    mem:"--mem=30G"
    tasks:"--ntasks=1"
    time:"24:00:00"
  }

  output {
    File trim_seq = "~{sequence_name}_trim.fastq.gz"
  }
}
