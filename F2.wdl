task create_alt_genome {
    File ref_genome

    command {
        /pirs/src/pirs/pirs diploid ${ref_genome} -s 0.001 -d 0 -v 0 -o alt
    }

    runtime {
        docker:"pirs-ddrad-cutadapt:v1"
    }

    output {
        File alt_fasta = "alt.snp.fa"
        File snps = "alt.snp.lst"
    }
}

task pedsim_files {
    File snp_file
    String genome_size
    String cmBymb

    command <<<
        R --vanilla --no-save <<RSCRIPT

            snps <- read.table("${snp_file}", stringsAsFactors = FALSE)
            n.marker <- dim(snps)[1]
            ## Map file
            # Marker names
            marker1 <- "M"
            marker2 <- 1:n.marker
            marker2 <- sprintf("%03d", marker2)
            marker <-paste0(marker1,marker2)
            # Chromossome and position
            tot = as.numeric("${genome_size}")*as.numeric("${cmBymb}")
            # The markers will be equally distribuited. There will be one marker each
            by = (tot)/(n.marker-1)
            pos <- seq(from = 0, to = tot, by = by)
            chr <- rep("C1",length(pos))
            map_file <- data.frame(marker=marker, chromosome=chr, position= pos)
            write.table(map_file, file = paste0("mapfile.map"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
            # Using sn1ps simulated by pirs genome described in ref_oryza.snp.lst file. Looking markers at the position:
            chr <- snps["V1"]
            pos <- snps["V2"]
            founder_file <- data.frame(marker=marker, P1_1=snps[["V4"]] , P1_2=snps[["V4"]], P2_1=snps[["V5"]], P2_2=snps[["V5"]])
            write.table(founder_file, file = paste0("founderfile.gen"), quote=FALSE, col.names = TRUE, row.names = FALSE, sep = "\t" )

            ## Parameters file
            parameter <- paste0("PLOIDY = 2
                                MAPFUNCTION = HALDANE
                                MISSING = NA
                                CHROMFILE = inb.chrom
                                POPTYPE = F2
                                POPSIZE = 150
                                MAPFILE = mapfile.map
                                FOUNDERFILE = founderfile.gen
                                OUTPUT = sim_inb")

            write.table(parameter, file = paste0("sim.par"), quote=FALSE, col.names = FALSE, row.names = FALSE, sep = "\t" )
            chrom <- data.frame("chromosome"= "C1", "length"= tot, "centromere"=tot/2, "prefPairing"= 0.0, "quadrivalents"=0.0)
            write.table(chrom, file= "inb.chrom", quote = F, col.names = T, row.names = F, sep= "\t")
        RSCRIPT
    >>>

    runtime {
        docker:"onemap:v1"
    }
    output {
        File mapfile="mapfile.map"
        File founderfile = "founderfile.gen"
        File parfile = "sim.par"
        File chromfile = "inb.chrom"
    }
}

task pedigreeSim{
    File map_file
    File founder_file
    File chrom_file
    File par_file

    command {
        sed -i 's+inb.chrom+${chrom_file}+g' ${par_file}
        sed -i 's+mapfile.map+${map_file}+g' ${par_file}
        sed -i 's+founderfile.gen+${founder_file}+g' ${par_file}

        java -jar /usr/jars/PedigreeSim.jar ${par_file}
    }

    runtime {
        docker:"java-in-the-cloud:v1"
    }

    output {
        File genotypes_dat = "sim_inb_genotypes.dat"
    }
}

task pedsim2vcf{
    File genotypes_dat
    File map_file
    File chrom_file
    File snp_file

    command <<<
        R --vanilla --no-save <<RSCRIPT
        
        library(onemap)              
        snps <- read.table("${snp_file}", stringsAsFactors = FALSE)
        pos <- snps[,3]
        chr <- snps[,1]

        pedsim2vcf(inputfile = "${genotypes_dat}", 
             map.file = "${map_file}", 
             chrom.file = "${chrom_file}",
             out.file = "simu.vcf",
             miss.perc = 0, counts = FALSE,pos = pos, haplo.ref = "P1_1", 
             chr = chr, phase = TRUE)

        RSCRIPT
    >>>
        
    runtime{
        docker: "onemap:v1"
    }

    output{
        File simu_vcf = "simu.vcf"
    }

}

# task sample_names{
#      File simu_vcf
#      command{
#      grep -i "CHROM" simu.vcf | cut -f1,2,3,4,5,6,7,8,9 --complement > id_names
#      tr -s '\t '  '\n'< id_names > id_names_lines
#      }
#      output{
#          File sampleNames = "id_names_lines"
#      }
#  }

task vcf2diploid{
    String sampleName
    File ref_genome
    File simu_vcf

    command{
        java -jar /usr/jars/vcf2diploid.jar -id ${sampleName} -chr ${ref_genome} -vcf ${simu_vcf}
    }
    runtime{
        docker:"java-in-the-cloud:v1"
    }
    output{
         File maternal_genomes = "Chr10_${sampleName}_maternal.fa"
         File paternal_genomes = "Chr10_${sampleName}_paternal.fa"
    }
}

task create_frags{
    String enzyme
    String sampleName
    File maternal_genomes
    File paternal_genomes

    command{
        /ddRADseqTools/Package/rsitesearch.py \
        --genfile=${maternal_genomes} \
        --fragsfile=${sampleName}_maternal_fragments.fasta \
        --rsfile=/ddRADseqTools/Package/restrictionsites.txt \
        --enzyme1=${enzyme} \
        --enzyme2=${enzyme} \
        --minfragsize=202 \
        --maxfragsize=500 \
        --fragstfile=${sampleName}_maternal_statistics.txt \
        --fragstinterval=25 \
        --plot=NO \
        --verbose=YES \
        --trace=NO

        /ddRADseqTools/Package/rsitesearch.py \
        --genfile=${paternal_genomes} \
        --fragsfile=${sampleName}_paternal_fragments.fasta \
        --rsfile=/ddRADseqTools/Package/restrictionsites.txt \
        --enzyme1=${enzyme} \
        --enzyme2=${enzyme} \
        --minfragsize=202 \
        --maxfragsize=500 \
        --fragstfile=${sampleName}_paternal_statistics.txt \
        --fragstinterval=25 \
        --plot=NO \
        --verbose=YES \
        --trace=NO

        cutadapt -l 202 \
        -o ${sampleName}_maternal_trim.fa \
        ${sampleName}_maternal_fragments.fasta

        cutadapt -l 202 \
        -o ${sampleName}_paternal_trim.fa \
        ${sampleName}_paternal_fragments.fasta
    }
    runtime{
        docker:"pirs-ddrad-cutadapt:v1"
    }
    output{
        File maternal_frags = "${sampleName}_maternal_fragments.fasta"
        File paternal_frags = "${sampleName}_paternal_fragments.fasta"
        File maternal_stats = "${sampleName}_maternal_statistics.txt"
        File paternal_stats = "${sampleName}_paternal_statistics.txt"
        File maternal_trim = "${sampleName}_maternal_trim.fa"
        File paternal_trim = "${sampleName}_paternal_trim.fa"
    }
}

task reads_simulations{
    File maternal_trim
    File paternal_trim
    String sampleName

    command{
        /pirs/src/pirs/pirs simulate --diploid ${maternal_trim} ${paternal_trim} \
        -l 100 -x 50 -m 150 -o ${sampleName}

        /cleanFastq/fixFastq "${sampleName}_100_150_1.fq" "${sampleName}_100_150_1_fix.fq"
        /cleanFastq/fixFastq "${sampleName}_100_150_2.fq" "${sampleName}_100_150_2_fix.fq"
    }
    runtime{
        docker:"pirs-ddrad-cutadapt:v1"
    }
    output{
        File reads1 = "${sampleName}_100_150_1_fix.fq"
        File reads2 = "${sampleName}_100_150_2_fix.fq"
    }
}

task alignment {
    command{
        bwa mem ${ref} ${reads1} ${reads2} > ${sampleName}.sam \
        samtools view -bS ${sampleName}.sam > ${sampleName}.bam \
        samtools sort ${sampleName}.bam -o ${sampleName}.sorted.bam \
        samtools index ${sampleName}.sorted.bam
    }
    runtime{
        docker:"bwa-samtools:v1"
    }
    output{
        input:
        File bam_file = "${sampleName}.sorted.bam"
    }
}

workflow F2 {
    
    File ref
    String genome_size
    String cmBymb
    File sampleNamesFile
    Array[String] sampleNames = read_lines(sampleNamesFile)
    String enzyme
    
    call create_alt_genome {
        input:
            ref_genome=ref 
    }

    call pedsim_files {
        input:
            snp_file = create_alt_genome.snps,
            genome_size=genome_size,
            cmBymb=cmBymb
    }

    call pedigreeSim {
        input:
            map_file = pedsim_files.mapfile,
            founder_file= pedsim_files.founderfile, 
            par_file=pedsim_files.parfile,
            chrom_file=pedsim_files.chromfile
    }

    call pedsim2vcf{
        input:
            map_file = pedsim_files.mapfile,
            chrom_file = pedsim_files.chromfile,
            genotypes_dat = pedigreeSim.genotypes_dat,
            snp_file = create_alt_genome.snps
    }

    scatter (sampleName in sampleNames) {
        call vcf2diploid{
            input: 
                sampleName = sampleName,
                ref_genome = ref,
                simu_vcf= pedsim2vcf.simu_vcf
        }

        call create_frags{
            input:
                sampleName = sampleName,
                enzyme = enzyme,
                maternal_genomes= vcf2diploid.maternal_genomes,
                paternal_genomes= vcf2diploid.paternal_genomes
        }

        call reads_simulations{
            input:
            sampleName = sampleName,
            maternal_trim = create_frags.maternal_trim,
            paternal_trim = create_frags.paternal_trim
        }
    }
}
