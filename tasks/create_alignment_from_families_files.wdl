version 1.0

import "alignment.wdl" as alg



struct Sample {
    String name
    Array[File] reads
    Array[String] libraries
}


struct Data {
    String experiment_name
    Array[String] names
    Array[Sample] samples
}

workflow CreateAlignmentFromFamilies {
    input {
        File families_info
        Reference references
        Int max_cores
    }

    call SepareIndividuals {
        input:
            families_info=families_info
    }

    scatter (sample in SepareIndividuals.dataset.samples) {
        call alg.RunBwaAlignment {
            input:
                sampleName  = sample.name,
                reads1      = sample.reads,
                libraries   = sample.libraries,
                references  = references,
                max_cores   = max_cores
        }
    }

    output {
        Array[Alignment] alignments = RunBwaAlignment.algn
        Array[File] bam = RunBwaAlignment.bam
        Array[File] bai = RunBwaAlignment.bai
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
        libs = {}

        with open("~{families_info}") as f:
            for l in f:
                file_path, sample_name, library = l.strip().split()
                sets.setdefault(sample_name, []).append(file_path)
                libs.setdefault(sample_name, []).append(library)

        names = [i for i in sets]
        samples = [{"name": sample, "reads": sets[sample], "libraries": libs[sample] } for sample in names]
        experiment_name = "Teste"
        print(json.dumps({"experiment_name": experiment_name, "samples":samples, "names": names}))
        CODE

    >>>

    runtime {
        docker: "python:3.7"
	mem:"--nodes=1"
	cpu:1
	time:"24:00:00"
    }

    output {
        Data dataset = read_json(stdout())
    }
}
