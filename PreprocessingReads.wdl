version 1.0

import "structs/preprocessing_reads_structs.wdl"

workflow PreprocessingReads{
    input {
      Specifications spec
    }

    Array[File] fq_files = read_lines(spec.raw_dict)

    call ProcessRadTags {
      input:
        enzyme = spec.enzyme,
        enzyme2 = spec.enzyme2,
        fq_files = fq_files,
        barcodes = spec.barcodes
    }

    scatter (sequence in ProcessRadTags.seq_results) {
      call RemoveAdapt {
        input:
          sequence = sequence,
          adapter = spec.adapter,
          sequence_name = basename(sequence)
      }
    }

    call TarFiles {
      input:
        sequences = RemoveAdapt.trim_seq
    }

    output {
      File results = TarFiles.results
    }
}


task ProcessRadTags {
    input {
      String enzyme
      String? enzyme2
      Array[File] fq_files
      File? barcodes
    }

    command <<<
      mkdir raw process_radtags_results
      mv ~{sep=" " fq_files} raw

      process_radtags -p raw/ -o process_radtags_results/ \
                      ~{"-b " + barcodes} \
                      --renz_1 ~{enzyme} ~{"--renz_2 " + enzyme2} \
                      -r -c -q -w 0.5

    >>>

    runtime {
      docker:"cristaniguti/stacks:0.0.1"
      job_name: "ProcessRadTags"
      node:"--nodes=1"
      mem:"--mem=30G"
      tasks:"--ntasks=1"
      time:"24:00:00"
    }

    output {
      Array[File] seq_results = glob("process_radtags_results/*.fq.gz")
    }
}


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

task TarFiles {
  input {
    Array[File] sequences
  }

  command <<<
    mkdir results
    mv ~{sep=" " sequences} results
    tar -czvf results.tar.gz results
  >>>

  runtime {
    docker:"kfdrc/cutadapt"
    job_name: "TarFiles"
    node:"--nodes=1"
    mem:"--mem=30G"
    tasks:"--ntasks=1"
    time:"24:00:00"
  }

  output {
    File results = "results.tar.gz"
  }
}
