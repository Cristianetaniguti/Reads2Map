version 1.0

import "structs/preprocessingS.wdl"

workflow PreprocessingReads{
    input{
      Specifications spec
    }

    call ProcessRadTags{
      input:
        enzyme = spec.enzyme,
        enzyme2 = spec.enzyme2,
        raw_dict = spec.raw_dict,
        barcodes = spec.barcodes
    }

    scatter(sequence in ProcessRadTags.seq_results){
      call RemoveAdapt{
        input:
          sequence = sequence,
          adapter = spec.adapter,
          sequence_name = basename(sequence)
      }
    }

    call TarFiles{
      input:
        sequences = RemoveAdapt.trim_seq
    }

    output{
      File results = TarFiles.results
    }
}


task ProcessRadTags{
    input{
      String enzyme
      String? enzyme2
      File raw_dict
      File? barcodes
    }

    command <<<
      mkdir raw process_radtags_results
      tar -xf ~{raw_dict} -C raw

      process_radtags -p raw/ -o process_radtags_results/ \
                      ~{"-b " + barcodes} \
                      --renz_1 ~{enzyme} ~{"--renz_2 " + enzyme2} \
                      -r -c -q -w 0.5

    >>>

    runtime{
      docker:"taniguti/stacks"
    }

    output{
      Array[File] seq_results = glob("process_radtags_results/*.fq.gz")
    }
}


task RemoveAdapt{
  input{
    File sequence
    String adapter
    String sequence_name
  }

  command <<<
    cutadapt -a ~{adapter} -o ~{sequence_name}_trim.fastq.gz ~{sequence} --minimum-length 64
  >>>

  runtime{
    docker:"kfdrc/cutadapt"
  }

  output{
    File trim_seq = "~{sequence_name}_trim.fastq.gz"
  }
}

task TarFiles{
  input{
    Array[File] sequences
  }

  command <<<
    mkdir results
    mv ~{sep=" " sequences} results
    tar -czvf results.tar.gz results
  >>>

  runtime{
    docker:"kfdrc/cutadapt"
  }

  output{
    File results = "results.tar.gz"
  }
}
