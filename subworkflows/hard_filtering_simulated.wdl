version 1.0

import "../structs/dna_seq_structs.wdl"
import "../tasks/gatk.wdl"
import "../tasks/utilsR.wdl"

workflow HardFiltering {
    input {
        ReferenceFasta references
        File vcf_file
        File vcf_tbi
        File? simu_vcf
        Int? seed
        Int? depth
    }

    call gatk.VariantsToTableForHardFilteringSimulated as VariantsToTable {
        input:
            vcf_file = vcf_file,
            vcf_tbi = vcf_tbi,
            simu_vcf = simu_vcf,
            reference = references.ref_fasta,
            reference_dict = references.ref_dict,
            reference_idx = references.ref_fasta_index
    }

    call utilsR.QualPlotsForHardFilteringSimulated as QualPlots {
        input:
            FalsePositives = VariantsToTable.FalsePositives,
            TruePositives  = VariantsToTable.TruePositives,
            Total          = VariantsToTable.Total,
            seed           = seed,
            depth          = depth
    }

    call gatk.VariantFiltration {
        input:
            vcf_file       = vcf_file,
            vcf_tbi        = vcf_tbi,
            reference      = references.ref_fasta,
            reference_idx  = references.ref_fasta_index,
            reference_dict = references.ref_dict
    }

    output {
        File Plots = QualPlots.Plots
        File filt_vcf  = VariantFiltration.filtered_vcf
    }
}
