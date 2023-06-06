version 1.0

import "../../structs/empirical_maps_structs.wdl"
import "../../structs/population_structs.wdl"

import "../../tasks/utils.wdl" as utils
import "../../tasks/utilsR.wdl" as utilsR
import "../../tasks/JointReports.wdl" as reports
import "../../tasks/mappoly.wdl" as mappoly_task

import "../../subworkflows/genotyping_empirical.wdl" as genotyping
import "../../subworkflows/snpcaller_maps_empirical.wdl" as snpcaller
import "../../subworkflows/gusmap_maps_empirical.wdl" as gusmap
import "../../subworkflows/mappoly_maps_empirical.wdl" as mappoly_sub

workflow Maps {

    input {
        Dataset dataset
        Array[File] vcfs
        Array[String] vcfs_software
        Array[String] vcfs_counts_source
        Boolean run_updog = true
        Boolean run_supermassa = false
        Boolean run_polyrad = true
        Boolean run_gusmap = false
        Boolean filter_noninfo
        String replaceADbyMissing 
        File? gatk_vcf_multi
        String gatk_mchap
        String? filters
        Int max_cores
        Int ploidy
        Float? prob_thres
        String? filt_segr
    }

    if (defined(filters)) {
        call utils.ApplyRandomFiltersArray {
            input:
                vcfs = vcfs,
                vcfs_SNPCall_software = vcfs_software,
                vcfs_Counts_source = vcfs_counts_source,
                vcfs_GenoCall_software = range(length(vcfs_software)),
                filters = filters
        }
    }

    Array[File] filtered_vcfs = select_first([ApplyRandomFiltersArray.vcfs_filt, vcfs])

    # Re-Genotyping with updog, supermassa and polyrad; and building maps with onemap
    scatter (idx in range(length(filtered_vcfs))) {

        call utils.SplitMarkers as splitgeno {
             input:
                vcf_file = filtered_vcfs[idx]
        }

        # Suggestion to improve performance of SuperMASSA, polyRAD and updog
        if(filter_noninfo){
            call utilsR.RemoveNonInformative {
            input:
                vcf_file = splitgeno.biallelics,
                parent1 = dataset.parent1,
                parent2 = dataset.parent2,
                replaceADbyMissing = replaceADbyMissing
            }
        }

        File vcf_up = select_first([RemoveNonInformative.vcf_filtered, splitgeno.biallelics])

        if(ploidy == 2) {
            if(run_updog){
                call genotyping.onemapMapsEmp as updogMaps {
                    input:
                        vcf_file = vcf_up,
                        SNPCall_program = vcfs_software[idx],
                        GenotypeCall_program = "updog",
                        CountsFrom = vcfs_counts_source[idx],
                        cross = dataset.cross,
                        parent1 = dataset.parent1,
                        parent2 = dataset.parent2,
                        chromosome = dataset.chromosome,
                        multiallelics = dataset.multiallelics,
                        multiallelics_file = splitgeno.multiallelics,
                        max_cores = max_cores,
                        ploidy = ploidy
                }
            }

            if(run_supermassa){
                call genotyping.onemapMapsEmp as supermassaMaps {
                    input:
                        vcf_file = vcf_up,
                        SNPCall_program = vcfs_software[idx],
                        GenotypeCall_program = "supermassa",
                        CountsFrom = vcfs_counts_source[idx],
                        cross = dataset.cross,
                        parent1 = dataset.parent1,
                        parent2 = dataset.parent2,
                        chromosome = dataset.chromosome,
                        multiallelics = dataset.multiallelics,
                        multiallelics_file = splitgeno.multiallelics,
                        max_cores = max_cores,
                        ploidy = ploidy
                }
            }

            if(run_polyrad){
                call genotyping.onemapMapsEmp as polyradMaps {
                    input:
                        vcf_file = vcf_up,
                        SNPCall_program = vcfs_software[idx],
                        GenotypeCall_program = "polyrad",
                        CountsFrom = vcfs_counts_source[idx],
                        cross = dataset.cross,
                        parent1 = dataset.parent1,
                        parent2 = dataset.parent2,
                        chromosome = dataset.chromosome,
                        multiallelics = dataset.multiallelics,
                        multiallelics_file = splitgeno.multiallelics,
                        max_cores = max_cores,
                        ploidy = ploidy
                }
            }
            
            if(run_gusmap){
                call gusmap.gusmapMapsEmp {
                        input:
                        vcf_file = vcf_up,
                        SNPCall_program = vcfs_software[idx],
                        CountsFrom = vcfs_counts_source[idx],
                        GenotypeCall_program = "gusmap",
                        parent1 = dataset.parent1,
                        parent2 = dataset.parent2,
                        max_cores = max_cores
                }
            }

            if(vcfs_counts_source[idx] != "bam"){
                call snpcaller.SNPCallerMapsEmp {
                    input:
                        vcf_file = splitgeno.biallelics,
                        cross = dataset.cross,
                        SNPCall_program = vcfs_software[idx],
                        GenotypeCall_program = "SNPCaller",
                        CountsFrom = vcfs_counts_source[idx],
                        parent1 = dataset.parent1,
                        parent2 = dataset.parent2,
                        chromosome = dataset.chromosome,
                        multiallelics = dataset.multiallelics,
                        max_cores = max_cores,
                        multiallelics_file = splitgeno.multiallelics,
                        multiallelics_mchap = gatk_vcf_multi,
                        mchap = gatk_mchap
                }
            }
        }

        if(ploidy > 2){
            if(run_updog){
                call mappoly_sub.MappolyMapsEmp as updogPolyMaps {
                    input:
                        vcf_file = vcf_up,
                        SNPCall_program = vcfs_software[idx],
                        GenotypeCall_program = "updog",
                        CountsFrom = vcfs_counts_source[idx],
                        cross = "F1",
                        parent1 = dataset.parent1,
                        parent2 = dataset.parent2,
                        max_cores = max_cores,
                        ploidy = ploidy,
                        prob_thres = prob_thres,
                        filt_segr = filt_segr
                }
            }

            if(run_polyrad){
                call mappoly_sub.MappolyMapsEmp as polyradPolyMaps {
                    input:
                        vcf_file = vcf_up,
                        SNPCall_program = vcfs_software[idx],
                        GenotypeCall_program = "polyrad",
                        CountsFrom = vcfs_counts_source[idx],
                        cross = "F1",
                        parent1 = dataset.parent1,
                        parent2 = dataset.parent2,
                        max_cores = max_cores,
                        ploidy = ploidy,
                        prob_thres = prob_thres,
                        filt_segr = filt_segr
                }
            }
            
            if(run_supermassa){
                call mappoly_sub.MappolyMapsEmp as supermassaPolyMaps {
                    input:
                        vcf_file = vcf_up,
                        SNPCall_program = vcfs_software[idx],
                        GenotypeCall_program = "supermassa",
                        CountsFrom = vcfs_counts_source[idx],
                        cross = "F1",
                        parent1 = dataset.parent1,
                        parent2 = dataset.parent2,
                        max_cores = max_cores,
                        ploidy = ploidy,
                        prob_thres = prob_thres,
                        filt_segr = filt_segr
                }
            }

            if(vcfs_counts_source[idx] != "bam" && vcfs_software[idx] != "stacks" && vcfs_software[idx] != "tassel"){
                call mappoly_task.MappolyReport {
                    input:
                        vcf_file = vcf_up,
                        SNPCall_program = vcfs_software[idx],
                        GenotypeCall_program = "SNPCaller",
                        CountsFrom = vcfs_counts_source[idx],
                        parent1 = dataset.parent1,
                        parent2 = dataset.parent2,
                        max_cores = max_cores,
			            ploidy = ploidy,
                        prob_thres = prob_thres,
                        filt_segr = filt_segr
                }
            }
        }
    }

    if(defined(updogMaps.tar_gz_report)) {
        Array[File] snpcaller_results = select_all(SNPCallerMapsEmp.tar_gz_report) 
    }
    if(defined(updogMaps.tar_gz_report)) {
        Array[File] updog_results =  select_all(updogMaps.tar_gz_report) 
    }
    if(defined(supermassaMaps.tar_gz_report)){
        Array[File] supermassa_results =  select_all(supermassaMaps.tar_gz_report)
    }
    if(defined(polyradMaps.tar_gz_report)){
        Array[File] polyrad_results = select_all(polyradMaps.tar_gz_report) 
    }
    if(defined(gusmapMapsEmp.tar_gz_report)){
        Array[File] gusmap_results = select_all(gusmapMapsEmp.tar_gz_report)
    }
    if(defined(MappolyReport.results)){
        Array[File] snpcaller_poly_results = select_all(MappolyReport.results)
    }
    if(defined(updogPolyMaps.tar_gz_report)){
        Array[File] updog_poly_results = select_all(updogPolyMaps.tar_gz_report)
    }
    if(defined(polyradPolyMaps.tar_gz_report)){
        Array[File] polyrad_poly_results = select_all(polyradPolyMaps.tar_gz_report)
    }
    if(defined(supermassaPolyMaps.tar_gz_report)){
        Array[File] supermassa_poly_results = select_all(supermassaPolyMaps.tar_gz_report)
    }
    
    # Compress files
    call reports.JointAllReports {
        input:
            SNPCallerMapsEmp = snpcaller_results,
            updogMaps = updog_results,
            polyradMaps = polyrad_results,
            supermassaMaps = supermassa_results,
            gusmapMapsEmp = gusmap_results,
            SNPCallerPolyMapsEmp = snpcaller_poly_results,
            updogPolyMaps = updog_poly_results,
            polyradPolyMaps = polyrad_poly_results,
            supermassaPolyMaps = supermassa_poly_results,
            max_cores = max_cores,
            ploidy = ploidy
    }

    output {
        File EmpiricalReads_results = JointAllReports.EmpiricalReads_results
    }
}
