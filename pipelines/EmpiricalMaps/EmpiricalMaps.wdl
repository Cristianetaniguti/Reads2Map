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
        Boolean filter_noninfo
        String replaceADbyMissing
        File? gatk_vcf_multi
        String gatk_mchap
        String? filters
        Int max_cores
        Int ploidy
    }

    if (defined(filters)) {
        call utils.ApplyRandomFiltersArray {
            input:
                vcfs = vcfs,
                vcfs_software = vcfs_software,
                vcfs_counts_source = vcfs_counts_source,
                filters = filters,
                chromosome = dataset.chromosome
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
            
            # Build maps with GUSMap
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
                    ploidy = ploidy
            }

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
                    ploidy = ploidy
            }

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
                    ploidy = ploidy
            }

            if(vcfs_counts_source[idx] != "bam"){
                call mappoly_task.MappolyReport {
                    input:
                        vcf_file = vcf_up,
                        SNPCall_program = vcfs_software[idx],
                        GenotypeCall_program = "SNPCaller",
                        CountsFrom = vcfs_counts_source[idx],
                        parent1 = dataset.parent1,
                        parent2 = dataset.parent2,
                        max_cores = max_cores,
			            ploidy = ploidy
                }
            }
        }
    }

    if(ploidy == 2){
        Array[File] snpcaller_results = select_all(SNPCallerMapsEmp.tar_gz_report)

        # Compress files
        call reports.JointReports {
            input:
                SNPCaller = snpcaller_results,
                updog = updogMaps.tar_gz_report,
                polyrad = polyradMaps.tar_gz_report,
                supermassa = supermassaMaps.tar_gz_report,
                gusmap = gusmapMapsEmp.tar_gz_report,
                max_cores = max_cores
        }
    }

    if(ploidy > 2){
        Array[File] snpcaller_results_poly = select_all(MappolyReport.results)

        call reports.JointReportsPoly {
            input:
                SNPCaller = snpcaller_results_poly,
                updog = updogPolyMaps.tar_gz_report,
                polyrad = polyradPolyMaps.tar_gz_report,
                supermassa = supermassaPolyMaps.tar_gz_report
        }
    }

    File Empirical_results_sele = select_first([JointReports.EmpiricalReads_results, JointReportsPoly.EmpiricalReads_results])

    output {
        File EmpiricalReads_results = Empirical_results_sele
    }
}
