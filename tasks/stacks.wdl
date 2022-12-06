version 1.0

task ProcessRadTags {
    input {
      String enzyme
      String? enzyme2
      Array[File] fq_files
      File? barcodes
    }

    Int disk_size = ceil(size(fq_files, "GiB") * 2)
    Int memory_size = 6000

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
      cpu:1
      # Cloud
      memory:"~{memory_size} MiB"
      disks:"local-disk " + disk_size + " HDD"
      # Slurm
      job_name: "ProcessRadTags"
      mem:"~{memory_size}M"
      time:"10:00:00"
    }

    output {
      Array[File] seq_results = glob("process_radtags_results/*.fq.gz")
    }
}
