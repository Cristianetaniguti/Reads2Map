version 1.0

import "./structs/reads_simuS.wdl"
import reads_simu.wdl as sub

workflow main{

    input{
        # não sei o que colocar aqui
    }

    Array[Int] seeds = read_lines("seeds_file") # pensei em adicionar um arquivo com os números dos seeds

    scatter(seed in seeds){
        call sub.reads_simu{
            input: 
            ReferenceFasta references # isso aqui esta errado
            Family family             # isso também
	        seed = seed
        }
    }

    call GraphicsAll{
        # vou fazer uma task para gerar gráficos dos resultados de todas as simulações
        # os inputs vao ser (só que de todos os seeds):
        input:
            cmBymb                    = family.cmBymb,
            depth                     = family.depth,
            mapfile                   = CreatePedigreeSimulatorInputs.mapfile,
            tot_mks                   = CreatePedigreeSimulatorInputs.tot_mks,
            gatkVCF_F                 = VcftoolsApplyFilters.gatkVCF_F,
            freebayesVCF_F            = VcftoolsApplyFilters.freebayesVCF_F,
            gatk_aval_vcf             = CalculateVcfMetrics.gatk_aval_vcf,
            freebayes_aval_vcf        = CalculateVcfMetrics.freebayes_aval_vcf,
            gatk_ref_depth            = CalculateVcfMetrics.gatk_ref_depth,
            gatk_ref_depth_bam        = BamCounts4Onemap.gatk_ref_bam,
            gatk_alt_depth            = CalculateVcfMetrics.gatk_alt_depth,
            gatk_alt_depth_bam        = BamCounts4Onemap.gatk_alt_bam,
            freebayes_ref_depth_bam   = BamCounts4Onemap.freebayes_ref_bam,
            freebayes_alt_depth_bam   = BamCounts4Onemap.freebayes_alt_bam,
            freebayes_ref_depth       = CalculateVcfMetrics.freebayes_ref_depth,
            freebayes_alt_depth       = CalculateVcfMetrics.freebayes_alt_depth,
            map_df                    = all_maps.map_df,
            map_GQ                    = all_maps.map_GQ,
            map_polyrad               = all_maps.map_polyrad,
            map_supermassa            = all_maps.map_supermassa,
            map_updog                 = all_maps.map_updog,
            map_bam_polyrad           = all_maps.map_bam_polyrad,
            map_bam_supermassa        = all_maps.map_bam_supermassa,
            map_bam_updog             = all_maps.map_bam_updog,
            error_info_GQ             = all_maps.error_info_GQ,
            error_info_updog          = all_maps.error_info_updog,
            error_info_polyrad        = all_maps.error_info_polyrad,
            error_info_supermassa     = all_maps.error_info_supermassa,
            error_info_bam_updog      = all_maps.error_info_updog,
            error_info_bam_polyrad    = all_maps.error_info_polyrad,
            error_info_bam_supermassa = all_maps.error_info_bam_supermassa
    }

    output{
        # Vão ser o gráficos e tabelas geradas pelo GraphicsAll
    }
}


task GraphicsAll{
    input{

    }

    command <<<
    >>>

    runtime{
        
    }

    output{

    }
}