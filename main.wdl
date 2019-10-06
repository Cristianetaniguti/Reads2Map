version 1.0

import "./structs/reads_simuS.wdl"
import "./tasks/reads_simu.wdl" as sub

workflow main{

    input{
        ReferenceFasta references
        FamilyTemplate family_template
        String name
        Int number_of_families
    }

    # ProduceFamiliesSeeds just generates random seeds. It returns an
    # array of integers
    call ProduceFamiliesSeeds {
        input:
            number_of_families=number_of_families
    }

    # Here we generate Family objects on the fly, based on the values
    # from the family_template and the random seed of the previous task.
    scatter(seed in ProduceFamiliesSeeds.seeds) {
        Family fam =  {
            "name": name,
            "cmBymb": family_template.cmBymb,
            "popsize": family_template.popsize,
            "enzyme": family_template.enzyme,
            "seed": seed,
            "depth": family_template.depth,
            "doses": family_template.doses,
            "ploidy": family_template.ploidy,
            "cross": family_template.cross
        }

        # Calling reads_simu for each seed
        call sub.reads_simu as ReadSimulations{
            input:
                references=references,
                family=fam
        }
    }

    # call GraphicsAll{
    #     # vou fazer uma task para gerar gráficos dos resultados de todas as simulações
    #     # os inputs vao ser (só que de todos os seeds):
    #     input:
    #         cmBymb                    = family.cmBymb,
    #         depth                     = family.depth,
    #         mapfile                   = CreatePedigreeSimulatorInputs.mapfile,
    #         tot_mks                   = CreatePedigreeSimulatorInputs.tot_mks,
    #         gatkVCF_F                 = VcftoolsApplyFilters.gatkVCF_F,
    #         freebayesVCF_F            = VcftoolsApplyFilters.freebayesVCF_F,
    #         gatk_aval_vcf             = CalculateVcfMetrics.gatk_aval_vcf,
    #         freebayes_aval_vcf        = CalculateVcfMetrics.freebayes_aval_vcf,
    #         gatk_ref_depth            = CalculateVcfMetrics.gatk_ref_depth,
    #         gatk_ref_depth_bam        = BamCounts4Onemap.gatk_ref_bam,
    #         gatk_alt_depth            = CalculateVcfMetrics.gatk_alt_depth,
    #         gatk_alt_depth_bam        = BamCounts4Onemap.gatk_alt_bam,
    #         freebayes_ref_depth_bam   = BamCounts4Onemap.freebayes_ref_bam,
    #         freebayes_alt_depth_bam   = BamCounts4Onemap.freebayes_alt_bam,
    #         freebayes_ref_depth       = CalculateVcfMetrics.freebayes_ref_depth,
    #         freebayes_alt_depth       = CalculateVcfMetrics.freebayes_alt_depth,
    #         map_df                    = all_maps.map_df,
    #         map_GQ                    = all_maps.map_GQ,
    #         map_polyrad               = all_maps.map_polyrad,
    #         map_supermassa            = all_maps.map_supermassa,
    #         map_updog                 = all_maps.map_updog,
    #         map_bam_polyrad           = all_maps.map_bam_polyrad,
    #         map_bam_supermassa        = all_maps.map_bam_supermassa,
    #         map_bam_updog             = all_maps.map_bam_updog,
    #         error_info_GQ             = all_maps.error_info_GQ,
    #         error_info_updog          = all_maps.error_info_updog,
    #         error_info_polyrad        = all_maps.error_info_polyrad,
    #         error_info_supermassa     = all_maps.error_info_supermassa,
    #         error_info_bam_updog      = all_maps.error_info_updog,
    #         error_info_bam_polyrad    = all_maps.error_info_polyrad,
    #         error_info_bam_supermassa = all_maps.error_info_bam_supermassa
    # }

    # Here you can reference outputs from the sub workflow. Remember that
    # it will be an array of the same type of the original.
    output {
        Array[File] coverages = ReadSimulations.coverage
    }
}


# task GraphicsAll{
#     input{

#     }

#     command <<<
#     >>>

#     runtime{

#     }

#     output{

#     }
# }

task ProduceFamiliesSeeds {
    input {
        Int number_of_families
    }

    command <<<
        python <<CODE
        import random
        for x in range(10):
            print(random.randint(1,101))
        CODE
    >>>

    runtime {
        docker: "python:3.7"
    }

    output {
        Array[Int] seeds = read_lines(stdout())
    }
}
