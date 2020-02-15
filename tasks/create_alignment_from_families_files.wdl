version 1.0

import "alignment.wdl" as alg
import "../structs/reads_simuS.wdl"



struct Sample {
    String name
    Array[File] reads
}


struct Data {
    String experiment_name
    Array[String] names
    Array[Sample] samples
}

workflow CreateAlignmentFromFamilies {
    input {
        File families_info
        ReferenceFasta references
    }

    call SepareIndividuals {
        input:
            families_info=families_info
    }

    scatter (sample in SepareIndividuals.dataset.samples) {
        call alg.RunBwaAlignment {
            input:
                sampleName = sample.name,
                reads1     = sample.reads,
                ref        = references.ref_fasta,
                geno_amb   = references.ref_amb,
                geno_ann   = references.ref_ann,
                geno_bwt   = references.ref_bwt,
                geno_pac   = references.ref_pac,
                geno_sa    = references.ref_sa
        }

        call alg.AddAlignmentHeader {
            input:
                sampleName = sample.name,
                bam_file   = RunBwaAlignment.bam_file,
                bam_idx    = RunBwaAlignment.bam_idx
        }
    }

    output {
        Array[Alignment] alignments = AddAlignmentHeader.algn
        Array[String] names = SepareIndividuals.dataset.names
    }
}


task SepareIndividuals {
    input {
        File families_info
    }

    command <<<
        python <<CODE
        import json
        sets = {}

        with open("~{families_info}") as f:
            for l in f:
                file_path, sample_name, library = l.strip().split()
                sets.setdefault(sample_name, []).append(file_path)

        samples = [{"name": k, "reads": v} for k, v in sets.items()]
        experiment_name = "Teste"
        print(json.dumps({"experiment_name": experiment_name, "samples":samples, "names": [i for i in sets]}))
        CODE
    >>>

    runtime {
        docker: "python:3.7"
    }

    output {
        Data dataset = read_json(stdout())
    }
}
