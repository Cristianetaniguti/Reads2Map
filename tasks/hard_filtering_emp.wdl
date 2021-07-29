version 1.0

import "../structs/snpcalling_empS.wdl"
import "../structs/reference_struct.wdl"
import "./split_filt_vcf.wdl" as norm_filt
import "./utils.wdl" as utils

workflow HardFilteringEmp {
    input {
        Reference references
        File vcf_file
        File vcf_tbi
    }

    call VariantsToTable {
        input:
            vcf_file = vcf_file,
            vcf_tbi = vcf_tbi,
            reference = references.ref_fasta,
            reference_dict = references.ref_dict,
            reference_idx = references.ref_fasta_index
    }

    call QualPlots {
        input:
            Total          = VariantsToTable.Total
    }
    
    call VariantFiltration {
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

task VariantsToTable {
    input {
        File vcf_file
        File vcf_tbi
        File reference
        File reference_dict
        File reference_idx
    }

    Int disk_size = ceil(size(reference, "GB") + size(vcf_file, "GB") + 2)

    command <<<

        gatk VariantsToTable \
            -V ~{vcf_file} \
            -F CHROM -F POS \
            -F QD \
            -F FS \
            -F SOR \
            -F MQ \
            -F MQRankSum \
            -F ReadPosRankSum \
            -O Total.table

    >>>

    runtime {
        docker: "taniguti/gatk-picard"
        # memory: "1 GB"
        # cpu: 1
        # disks: "local-disk " + disk_size + " HDD"
        job_name: "VariantsToTable"
        node:"--nodes=1"
        mem:"--mem=10G"
        tasks:"--ntasks=1"
        time:"00:50:00"
    }

    output {
        File Total = "Total.table"
    }

}


task QualPlots {
    input {
        File Total
    }

    Int disk_size = ceil(size(Total, "GB") + 1)

    command <<<
        R --vanilla --no-save <<RSCRIPT
            system("cp ~{Total} .")

            library(ggplot2)
            library(dplyr)
            library(tidyr)

            tot <- read.table("Total.table", header = T)
            tot <- cbind(set = "Total", tot)

            df <- pivot_longer(tot, cols = c(4:9))

            p <- df %>% filter(name == "QD") %>%   
                    ggplot(aes(x=value, fill=set)) +
                    geom_density(alpha=0.4) +
                    geom_vline(aes(xintercept=2), color = "purple", linetype="dashed") +
                    xlab("QD")

            ggsave(p, filename = "QD.png")

            p <- df %>% filter(name == "FS") %>%   
                    ggplot(aes(x=value, fill=set)) +
                        geom_density(alpha=0.4) +
                        xlab("FS") + geom_vline(aes(xintercept=60), color = "purple", linetype="dashed") 

            ggsave(p, filename = "FS.png")

            p <- df %>% filter(name == "SOR") %>%   
                    ggplot(aes(x=value, fill=set)) +
                        geom_density(alpha=0.4) +
                        xlab("SOR") +
                        geom_vline(aes(xintercept=3), color = "purple", linetype="dashed") 

            ggsave(p, filename = "SOR.png")

            p <- df %>% filter(name == "MQ") %>%   
                    ggplot(aes(x=value, fill=set)) +
                        geom_density(alpha=0.4) +
                        xlab("MQ") +
                        geom_vline(aes(xintercept=40), color = "purple", linetype="dashed") 

            ggsave(p, filename = "MQ.png")

            p <- df %>% filter(name == "MQRankSum") %>%   
                    ggplot(aes(x=value, fill=set)) +
                        geom_density(alpha=0.4) +
                        xlab("MQRankSum") +
                        geom_vline(aes(xintercept=-12.5), color = "purple", linetype="dashed") 

            ggsave(p, filename = "MQRankSum.png")

            p <- df %>% filter(name == "ReadPosRankSum") %>%   
                    ggplot(aes(x=value, fill=set)) +
                        geom_density(alpha=0.4) +
                        xlab("ReadPosRankSum") +
                        geom_vline(aes(xintercept=-8), color = "purple", linetype="dashed") 

            ggsave(p, filename = "ReadPosRankSum.png")
        
            system("mkdir QualPlots")
            system("mv *png QualPlots")
            system("tar -czvf QualPlots.tar.gz QualPlots")

        RSCRIPT
        
    >>>

    runtime {
        docker: "cristaniguti/reads2map"
        # memory: "1 GB"
        # cpu: 1
        # disks: "local-disk " + disk_size + " HDD"
        job_name: "QualPlots"
        node:"--nodes=1"
        mem:"--mem=10G"
        tasks:"--ntasks=1"
        time:"00:20:00"
    }

    output {
        File Plots = "QualPlots.tar.gz"
    }

}


task VariantFiltration {
    input {
        File vcf_file
        File vcf_tbi
        File reference
        File reference_idx
        File reference_dict
    }
    Int disk_size = ceil(size(vcf_file, "GB") + size(reference, "GB") + 1)

    command <<<
        /gatk/gatk VariantFiltration \
            -V ~{vcf_file} \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O gatk_filters.vcf.gz

        /gatk/gatk SelectVariants \
            -R ~{reference} \
            -V gatk_filters.vcf.gz \
            --exclude-filtered \
            -O gatk_filtered.vcf.gz
            
    >>>

    runtime {
        docker: "taniguti/gatk-picard"
        # memory: "1 GB"
        # cpu: 1
        # disks: "local-disk " + disk_size + " HDD"
        job_name: "VariantFiltration"
        node:"--nodes=1"
        mem:"--mem=10G"
        tasks:"--ntasks=1"
        time:"00:10:00"
    }

    output {
        File filters_vcf = "gatk_filters.vcf.gz"
        File filtered_vcf = "gatk_filtered.vcf.gz"
    }
}
