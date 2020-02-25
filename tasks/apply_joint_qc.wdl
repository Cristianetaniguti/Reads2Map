version 1.0

import "../structs/alignment_struct.wdl"
import "./utils.wdl" as utils


workflow ApplyJointQC {
    input {
      String program
      Array[File] vcfs
      Array[File] tbis
    }

    call utils.VcftoolsMerge {
      input:
        prefix=program,
        vcfs=vcfs,
        tbis=tbis
    }

    call utils.VcftoolsApplyFilters {
      input:
        vcf_in=VcftoolsMerge.vcf,
        tbi_in=VcftoolsMerge.tbi,
        max_missing=0.75,
        min_alleles=2,
        max_alleles=2,
        maf=0.05,
        program=program
    }

    output {
      File vcf = VcftoolsApplyFilters.vcf
      File tbi = VcftoolsApplyFilters.tbi
    }
}
