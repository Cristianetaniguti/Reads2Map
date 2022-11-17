version 1.0


import "../../structs/preprocessing_reads_structs.wdl"

import "../../tasks/stacks.wdl"
import "../../tasks/cutadapt.wdl"
import "../../tasks/utils.wdl"


workflow PreprocessingReads{
    input {
      Specifications spec
    }

    Array[File] fq_files = read_lines(spec.raw_dict)

    call stacks.ProcessRadTags {
      input:
        enzyme = spec.enzyme,
        enzyme2 = spec.enzyme2,
        fq_files = fq_files,
        barcodes = spec.barcodes
    }

    scatter (sequence in ProcessRadTags.seq_results) {
      call cutadapt.RemoveAdapt {
        input:
          sequence = sequence,
          adapter = spec.adapter,
          sequence_name = basename(sequence)
      }
    }

    call utils.TarFiles {
      input:
        sequences = RemoveAdapt.trim_seq
    }

    output {
      File results = TarFiles.results
    }
}
