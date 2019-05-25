task create_alt_genome {
    File ref_genome

     command {
        /pirs/src/pirs/pirs diploid ${ref_genome} -s 0.001 -d 0 -v 0 -o alt
        chmod 777 "alt.snp.lst" "alt.snp.fa"
     }

     runtime {
        docker:"pirs"
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
        docker:"r-base:3.6.0"
    }
    output {
        File mapfile="mapfile.map"
        File founderfile = "founderfile.gen"
        File parfile = "sim.par"
        File chromfile = "inb.chrom"
    }
}

task pedigreeSim{
    File pedigreeSimJar
    File map_file
    File founder_file
    File chrom_file
    File par_file

    command {
        sed -i 's+inb.chrom+${chrom_file}+g' ${par_file}
        sed -i 's+mapfile.map+${map_file}+g' ${par_file}
        sed -i 's+founderfile.gen+${founder_file}+g' ${par_file}

        java -jar ${pedigreeSimJar} ${par_file} 
    }
    output {
        File genotypes_dat = "sim_inb_genotypes.dat"
    }
}

workflow F2 {
    
    File ref
    String genome_size
    String cmBymb
    File pedigreeSim_jar
    
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
            chrom_file=pedsim_files.chromfile,
            pedigreeSimJar=pedigreeSim_jar
    }
}