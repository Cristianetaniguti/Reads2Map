version 1.0

import "../../structs/empirical_maps_structs.wdl"
import "../../structs/population_structs.wdl"

import "../../tasks/utils.wdl" as utils
import "../../tasks/utilsR.wdl" as utilsR
import "../../tasks/JointReports.wdl" as reports

import "../../subworkflows/genotyping_empirical.wdl" as genotyping
import "../../subworkflows/snpcaller_maps_empirical.wdl" as snpcaller
import "../../subworkflows/gusmap_maps_empirical.wdl" as gusmap

workflow Maps {

    input {
        Dataset dataset
        File? gatk_vcf_multi
        String gatk_mchap
        File? gatk_vcf
        File? freebayes_vcf
        File? gatk_vcf_bam_counts
        File? freebayes_vcf_bam_counts
        String? filters
        Int max_cores
    }

    if (defined(filters)) {
        call utils.ApplyRandomFilters {
            input:
                gatk_vcf = gatk_vcf,
                freebayes_vcf = freebayes_vcf,
                gatk_vcf_bam_counts = gatk_vcf_bam_counts,
                freebayes_vcf_bam_counts = freebayes_vcf_bam_counts,
                filters = filters,
                chromosome = dataset.chromosome
        }
    }

    File filtered_gatk_vcf = select_first([ApplyRandomFilters.gatk_vcf_filt, gatk_vcf])
    File filtered_gatk_vcf_bamcounts = select_first([ApplyRandomFilters.gatk_vcf_bam_counts_filt, gatk_vcf_bam_counts])
    File filtered_freebayes_vcf = select_first([ApplyRandomFilters.freebayes_vcf_filt, freebayes_vcf])
    File filtered_freebayes_vcf_bamcounts = select_first([ApplyRandomFilters.freebayes_vcf_bam_counts_filt, freebayes_vcf_bam_counts])

    PopulationAnalysis gatk_processing = {"method": "gatk", "vcf": filtered_gatk_vcf, "bam": filtered_gatk_vcf_bamcounts}
    PopulationAnalysis freebayes_processing = {"method": "freebayes", "vcf": filtered_freebayes_vcf, "bam": filtered_freebayes_vcf_bamcounts}

    # Re-Genotyping with updog, supermassa and polyrad; and building maps with onemap
    scatter (analysis in [gatk_processing, freebayes_processing]) {

        Map[String, File] vcfs = {"vcf": analysis.vcf, "bam": analysis.bam}

        scatter (origin in ["vcf", "bam"]) {

            call utils.SplitMarkers as splitgeno {
                 input:
                    vcf_file = vcfs[origin]
            }

            # Suggestion to improve performance of SuperMASSA, polyRAD and updog
            call utilsR.FilterSegregation {
             input:
                vcf_file = splitgeno.biallelics,
                parent1 = dataset.parent1,
                parent2 = dataset.parent2
            }

            call genotyping.onemapMapsEmp as updogMaps {
                input:
                    vcf_file = FilterSegregation.vcf_filtered,
                    # vcf_file = splitgeno.biallelics,
                    SNPCall_program = analysis.method,
                    GenotypeCall_program = "updog",
                    CountsFrom = origin,
                    cross = dataset.cross,
                    parent1 = dataset.parent1,
                    parent2 = dataset.parent2,
                    chromosome = dataset.chromosome,
                    multiallelics = dataset.multiallelics,
                    multiallelics_file = splitgeno.multiallelics,
                    max_cores = max_cores
            }

            call genotyping.onemapMapsEmp as supermassaMaps {
                input:
                    vcf_file = FilterSegregation.vcf_filtered,
                    # vcf_file = splitgeno.biallelics,
                    SNPCall_program = analysis.method,
                    GenotypeCall_program = "supermassa",
                    CountsFrom = origin,
                    cross = dataset.cross,
                    parent1 = dataset.parent1,
                    parent2 = dataset.parent2,
                    chromosome = dataset.chromosome,
                    multiallelics = dataset.multiallelics,
                    multiallelics_file = splitgeno.multiallelics,
                    max_cores = max_cores
            }

            call genotyping.onemapMapsEmp as polyradMaps {
                input:
                    vcf_file = FilterSegregation.vcf_filtered,
                    # vcf_file = splitgeno.biallelics,
                    SNPCall_program = analysis.method,
                    GenotypeCall_program = "polyrad",
                    CountsFrom = origin,
                    cross = dataset.cross,
                    parent1 = dataset.parent1,
                    parent2 = dataset.parent2,
                    chromosome = dataset.chromosome,
                    multiallelics = dataset.multiallelics,
                    multiallelics_file = splitgeno.multiallelics,
                    max_cores = max_cores
            }
        }

       call utils.SplitMarkers as splitvcf {
            input:
              vcf_file = analysis.vcf
       }

       call utils.SplitMarkers as splitbam {
            input:
              vcf_file = analysis.bam
       }

       # Build maps with GUSMap
       call gusmap.gusmapMapsEmp {
            input:
              vcf_file = splitvcf.biallelics,
              vcf_bam_file = splitbam.biallelics,
              SNPCall_program = analysis.method,
              GenotypeCall_program = "gusmap",
              parent1 = dataset.parent1,
              parent2 = dataset.parent2,
              max_cores = max_cores
        }

        call snpcaller.SNPCallerMapsEmp {
            input:
                vcf_file = splitvcf.biallelics,
                cross = dataset.cross,
                SNPCall_program = analysis.method,
                GenotypeCall_program = "SNPCaller",
                CountsFrom = "vcf",
                parent1 = dataset.parent1,
                parent2 = dataset.parent2,
                chromosome = dataset.chromosome,
                multiallelics = dataset.multiallelics,
                max_cores = max_cores,
                multiallelics_file = splitvcf.multiallelics,
                multiallelics_mchap = gatk_vcf_multi,
                mchap = gatk_mchap
        }
    }

    # Compress files
    call reports.JointReports {
        input:
            SNPCaller = SNPCallerMapsEmp.tar_gz_report,
            updog = flatten(updogMaps.tar_gz_report),
            polyrad = flatten(polyradMaps.tar_gz_report),
            supermassa = flatten(supermassaMaps.tar_gz_report),
            gusmap = gusmapMapsEmp.tar_gz_report,
            max_cores = max_cores
    }

    output {
        File EmpiricalReads_results = JointReports.EmpiricalReads_results
    }
}
