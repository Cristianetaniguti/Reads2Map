version 1.0

import "structs/snpcalling_empS.wdl"
import "structs/maps_empS.wdl"

import "tasks/snpcalling_emp.wdl" as snpcalling
import "tasks/maps_emp.wdl" as maps

workflow EmpiricalReads {

    input {
        Samples_info samples_info
        Reference references
        Dataset dataset
        SplitVCF splitvcf
        String? filters
        Int max_cores
    }

    # TODO: Conferir splitvcf
    call snpcalling.SNPCalling {
        input:
            samples_info = samples_info,
            references = references,
            splitvcf = splitvcf,
            max_cores = max_cores
    }

    call maps.Maps {
        input:
            dataset = dataset,
            gatk_vcf = SNPCalling.gatk_vcf_bi,
            freebayes_vcf = SNPCalling.freebayes_vcf_bi,
            gatk_multi = SNPCalling.gatk_vcf_multi,
            freebayes_multi = SNPCalling.freebayes_vcf_multi,
            gatk_vcf_bam_counts = SNPCalling.gatk_vcf_bi_bam_count,
            freebayes_vcf_bam_counts = SNPCalling.freebayes_vcf_bi_bam_count,
            filters = filters
    }

    output {
        File EmpiricalReads_results = Maps.EmpiricalReads_results
        File gatk_vcf_bi = SNPCalling.gatk_vcf_bi
        File gatk_vcf_multi = SNPCalling.gatk_vcf_multi
        File freebayes_vcf_bi = SNPCalling.freebayes_vcf_bi
        File freebayes_vcf_multi = SNPCalling.freebayes_vcf_multi
    }
}
