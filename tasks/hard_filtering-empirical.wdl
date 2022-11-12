version 1.0

import "../structs/struct_reference.wdl"
import "gatk.wdl"
import "custom/r_libs.wdl"

workflow HardFilteringEmp {
    input {
        Reference references
        File vcf_file
        File vcf_tbi
    }

    call gatk.VariantsToTable {
        input:
            vcf_file = vcf_file,
            vcf_tbi = vcf_tbi,
            reference = references.ref_fasta,
            reference_dict = references.ref_dict,
            reference_idx = references.ref_fasta_index
    }

    call r_libs.QualPlots {
        input:
            Total=VariantsToTable.Total
    }

    call gatk.VariantFiltration {
        input:
            vcf_file = vcf_file,
            vcf_tbi = vcf_tbi,
            reference = references.ref_fasta,
            reference_idx  = references.ref_fasta_index,
            reference_dict = references.ref_dict
    }

    output {
        File Plots = QualPlots.Plots
        File filt_vcf = VariantFiltration.filtered_vcf
    }
}
