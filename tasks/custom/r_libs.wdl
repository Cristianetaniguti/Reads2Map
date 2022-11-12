version 1.0

import "../../structs/struct_reference.wdl"

task SimuscopProfile{
  input {
    String library_type
    File?  emp_bam
    File   vcf
    Reference  references
  }

  Int disk_size = ceil(size(vcf, "GiB") * 2 + size(emp_bam, "GiB")) 
  Int memory_size = 5000

  command <<<
    R --vanilla --no-save <<RSCRIPT

    library(vcfR)
    library(simuscopR)

    if("~{library_type}" == "exome"){

      system(paste0("bamtobed -i ~{emp_bam} > bed_file"))

      seqToProfile("~{emp_bam}", "bed_file", "~{vcf}",
             "~{references.ref_fasta}", "profile")

    } else {
      seqToProfile("~{emp_bam}", vcf.file =  "~{vcf}",
             reference = "~{references.ref_fasta}",  out.profile = "sample.profile")
    }

    RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/reads2map:0.0.1"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "SimuscopProfile"
    mem:"~{memory_size}M"
    time:"05:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Run [seqToProfile](https://github.com/qasimyu/simuscop) to generate data set profile."
  }

  output{
    File profile = "sample.profile"
  }
}