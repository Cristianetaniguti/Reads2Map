version 1.0

import "../tasks/stacks.wdl"


workflow StacksGenotyping {
    input {
        Array[File] bams
        File? pop_map
        Int max_cores
    }

    if(!defined(pop_map)){
        call stacks.CreatePopMapFile {
            input:
                bams = bams
        }
    }

    File pop_map_sele = select_first([pop_map, CreatePopMapFile.pop_map])
    
    call stacks.RefMap {
        input:
            bams = bams,
            pop_map = pop_map_sele,
            max_cores = max_cores
    }

    output {
        Array[File] vcfs = [RefMap.stacks_vcf]
        File stacks_multiallelics = RefMap.stacks_multiallelics
        Array[String] software_sele = ["stacks"]
        Array[String] source_sele = ["vcf"]
    }
}