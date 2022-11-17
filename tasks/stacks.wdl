version 1.0

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
