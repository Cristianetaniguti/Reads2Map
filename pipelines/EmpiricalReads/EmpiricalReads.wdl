version 1.0

import "../../structs/empirical_maps_structs.wdl"
import "../../structs/dna_seq_structs.wdl"

import "../../pipelines/EmpiricalSNPCalling/EmpiricalSNPCalling.wdl" as snpcalling
import "../../pipelines/EmpiricalMaps/EmpiricalMaps.wdl" as maps

workflow EmpiricalReads {

    input {
        File samples_info
        ReferenceFasta references
        Dataset dataset
        Int max_cores
        Int chunk_size
        Boolean rm_dupli = true
        Boolean gatk_mchap = false
        Boolean hardfilters = true
        Boolean replaceAD = true
        Boolean run_gatk = true
        Boolean run_freebayes = true
        Int ploidy = 2
        Int n_chrom
        String? filters
    }

    call snpcalling.SNPCalling {
        input:
            samples_info = samples_info,
            references = references,
            max_cores = max_cores,
            rm_dupli = rm_dupli,
            P1 = dataset.parent1,
            P2 = dataset.parent2,
            gatk_mchap = gatk_mchap,
            hardfilters = hardfilters,
            replaceAD = replaceAD,
            run_gatk = run_gatk,
            run_freebayes = run_freebayes,
            ploidy = ploidy,
            n_chrom = n_chrom
    }

    call maps.Maps {
        input:
            dataset = dataset,
            gatk_vcf_multi = SNPCalling.gatk_multi_vcf,
            gatk_mchap = gatk_mchap,
            gatk_vcf = SNPCalling.gatk_vcf,
            freebayes_vcf = SNPCalling.freebayes_vcf,
            gatk_vcf_bam_counts = SNPCalling.gatk_vcf_bam_count,
            freebayes_vcf_bam_counts = SNPCalling.freebayes_vcf_bam_count,
            filters = filters,
            max_cores = max_cores
    }

    output {
        File EmpiricalReads_results = Maps.EmpiricalReads_results
    }
}