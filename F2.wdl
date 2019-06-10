# Tenho genoma de referencia e vou criar o homologo dele
# TODO: simular indel
task create_alt_genome {
    File ref_genome

    command {
        # taxa de snp: 0.001
        /pirs/src/pirs/pirs diploid ${ref_genome} -s 0.001 -d 0 -v 0 -o alt --random-seed 1515
    }

    runtime {
        docker:"pirs-ddrad-cutadapt:v1"
    }

    output {
        File alt_fasta = "alt.snp.fa" # homologo com snps
        File snps = "alt.snp.lst"  # onde criou os snps
    }
}

# Descreve quem sao os pais, com os polimorfismos definidos (lista anterior)
# Output: simula uma progene (simula recombinacoes)
# Aqui simula o mapa de ligacao
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
                                SEED = 8899
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

# Apenas executa o programa com parametros acima
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

# Script para transformar .dat em VCF
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

# pega os polimorfismos do VCF e transfere para o genoma.
# Lida bem com indels tb
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
         File maternal_genomes = "Chr10_${sampleName}_maternal.fa" # Cada fasta representa um cromossomo
         File paternal_genomes = "Chr10_${sampleName}_paternal.fa"
    }
}

# Tecnica usada para obter mapas. Recorta o genoma de acordo com
# alguma enzima de restricao e sequencia X bases iniciais. Nesta task
# foi definido 202 pb.
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
        File maternal_frags = "${sampleName}_maternal_fragments.fasta"  # fragmentos obtidos ao cortar com a enzima
        File paternal_frags = "${sampleName}_paternal_fragments.fasta"
        File maternal_stats = "${sampleName}_maternal_statistics.txt"
        File paternal_stats = "${sampleName}_paternal_statistics.txt"
        File maternal_trim = "${sampleName}_maternal_trim.fa"  # apos o trim pelo cutadapt (202 pb)
        File paternal_trim = "${sampleName}_paternal_trim.fa"
    }
}

# Utiliza os reads trim da task anterior para simular o sequenciamento Illumina
# TODO: simular single end apenas
# -l:
# -x: profundidade do sequenciamento (devera ser parametrizavel)
# -m:
# fixFastq para acertar a qualidade. As vezes os ultimos reads
# não apresentavam as respectivas qualidades.
task reads_simulations{
    File maternal_trim
    File paternal_trim
    String sampleName

    command{
        /pirs/src/pirs/pirs simulate \
            --diploid ${maternal_trim} ${paternal_trim} \
            -l 100 -x 50 -m 150 -o ${sampleName} --random-seed 1515

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
        String sampleName
        File ref
        File reads1
        File reads2
        File geno_amb
        File geno_ann
        File geno_bwt
        File geno_pac
        File geno_sa

    command{
        /usr/gitc/bwa mem ${ref} ${reads1} ${reads2} | \
            java -jar /usr/gitc/picard.jar SortSam \
                I=/dev/stdin \
                O=${sampleName}.sorted.bam \
                SORT_ORDER=coordinate \
                CREATE_INDEX=true
        mv ${sampleName}.sorted.bai ${sampleName}.sorted.bam.bai
    }
    runtime{
        docker:"us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
    }
    output{
        File bam_file = "${sampleName}.sorted.bam"
	    File bam_idx = "${sampleName}.sorted.bam.bai"
    }
}

task add_labs{
    String sampleName
    File bam_file
    File bam_idx

    command{
        mkdir tmp
        java -jar /gatk/picard.jar AddOrReplaceReadGroups \
            I=${bam_file} \
            O=${sampleName}_rg.bam \
            RGLB=lib-${sampleName} \
            RGPL=illumina \
            RGID=FLOWCELL1.LANE1.${sampleName} \
            RGSM=${sampleName} \
            RGPU=FLOWCELL1.LANE1.${sampleName} \
            CREATE_INDEX=true \
            TMP_DIR=tmp

        mv ${sampleName}_rg.bai ${sampleName}_rg.bam.bai
    }
    runtime{
        docker:"gatk-picard:v1"
    }
    output{
        File bam_rg = "${sampleName}_rg.bam"
        File bam_rg_index = "${sampleName}_rg.bam.bai"
    }
}


task HaplotypeCallerERC {

    File ref
    File geno_fai
    String sampleName
    File bam_rg
    File bam_rg_idx
    File geno_dict

    command {
        /gatk/gatk HaplotypeCaller \
            -ERC GVCF \
            -R ${ref} \
            -I ${bam_rg} \
            -O ${sampleName}_rawLikelihoods.g.vcf
    }
    runtime{
        docker:"gatk-picard:v1"
    }
    output {
        File GVCF = "${sampleName}_rawLikelihoods.g.vcf"
        File GVCF_idx = "${sampleName}_rawLikelihoods.g.vcf.idx"
    }
}

task create_gatk_database {
    String path_gatkDatabase
    Array[File] GVCFs
    Array[File] GVCFs_idx
    
    command {
        /gatk/gatk GenomicsDBImport \
            --genomicsdb-workspace-path ${path_gatkDatabase} \
            -L Chr10 \
            -V ${sep=" -V " GVCFs} 

        tar -cf ${path_gatkDatabase}.tar ${path_gatkDatabase}
    }
    runtime {
        docker:"gatk-picard:v1"
    }
    output {
        File workspace_tar = "${path_gatkDatabase}.tar"
    }
}

task GenotypeGVCFs {
    
    File workspace_tar
    String output_vcf_filename
    File ref
    File geno_fai
    File geno_dict


    command{
        tar -xf ${workspace_tar}
        WORKSPACE=$( basename ${workspace_tar} .tar)

        /gatk/gatk  GenotypeGVCFs \
            -R ${ref} \
            -O ${output_vcf_filename} \
            -G StandardAnnotation \
            -V gendb://$WORKSPACE 
    }

    runtime {
        docker: "gatk-picard:v1"
    }
    output {
        File output_vcf = "${output_vcf_filename}"
        File output_vcf_index = "${output_vcf_filename}.idx"
    }
}

task freebayes {
    String freebayesVCFname
    File ref
    Array[File] bam_rg
    command{
        freebayes -f ${ref} ${sep=" " bam_rg} > ${freebayesVCFname}
    }
    runtime{
        docker:"freebayes:v1"
    }
    output{
        File freebayesVCF = "${freebayesVCFname}"
    }
}

# process_radtags only discards by absense of cut site
# task process_radtags{
#     String sampleName
#     File reads1
#     File reads2
#     String enzyme
#     command{
#         process_radtags -1 ${reads1} -2 ${reads2}  -o . -e ${enzyme} 
#     }
#     rumtime{
#         docker:"stacks:v1"
#     }
#     output{
#         File process_reads1 = "${sampleName}_100_150_1_fix.1.fq"
#         File process_reads1 = "${sampleName}_100_150_2_fix.2.fq"
#     }
# }

# apenas prepara arquivo
task create_popmapFile{
    File sampleNamesFile
    
    command <<<
        R --vanilla --no-save <<RSCRIPT
        names <- read.table("${sampleNamesFile}", header = F)
        mapnames <- paste0(t(names), '_rg')
        mapdf <- data.frame(mapnames, rep(1, length(mapnames)))
        write.table(mapdf, file = 'popmap.txt', sep = '\t', col.names = F, row.names = F, quote=F)

        RSCRIPT
    >>>

    runtime{
        docker:"onemap:v1"
    }
    output{
        File popmapfile = "popmap.txt"
    }
}

# snp calling (analogo a gatk e freebayes)
task ref_map {
    Array[File] bam_rg
    File popmapfile
    
    command{
        cp ${sep=" " bam_rg} .
        ref_map.pl --samples . --popmap ${popmapfile} -o . -X "populations:--vcf"
    }
    runtime{
        docker:"stacks:v1"
    }
    output{
        File stacksVCF = stdout()
    }
}


workflow F2 {
    
    String genome_size
    String cmBymb
    File sampleNamesFile
    Array[String] sampleNames = read_lines(sampleNamesFile)
    String enzyme
    String path_gatkDatabase
    String gatkVCFname
    String freebayesVCFname  
    String enzymeName

    File ref
    File ref_dict
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File ref_fai

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
	
	    call alignment {
	        input:
		        sampleName = sampleName,
                reads1 = reads_simulations.reads1,
                reads2 = reads_simulations.reads2,
                ref = ref,
                geno_amb = ref_amb,
     	        geno_ann = ref_ann,
                geno_bwt = ref_bwt,
                geno_pac = ref_pac,
                geno_sa	 = ref_sa	 
	    }

        call add_labs{
            input:
                sampleName = sampleName,
                bam_file=alignment.bam_file,
                bam_idx = alignment.bam_idx
        }

        call HaplotypeCallerERC{
            input:
                sampleName = sampleName,
                ref = ref,
                geno_fai = ref_fai,
                bam_rg = add_labs.bam_rg,
                bam_rg_idx = add_labs.bam_rg_index,
                geno_dict = ref_dict
        }
    }
    call create_gatk_database{
        input:
            GVCFs = HaplotypeCallerERC.GVCF,
            GVCFs_idx = HaplotypeCallerERC.GVCF_idx,
            path_gatkDatabase = path_gatkDatabase
    }

    call GenotypeGVCFs {
        input:
            workspace_tar=create_gatk_database.workspace_tar,
            output_vcf_filename = gatkVCFname,
            ref = ref,
            geno_fai = ref_fai,
            geno_dict = ref_dict
    }

    call freebayes{
         input:
            ref = ref,
            bam_rg = add_labs.bam_rg,
            freebayesVCFname = freebayesVCFname
    }
    
    call create_popmapFile{
        input:
            sampleNamesFile = sampleNamesFile
    }

    call ref_map{
        input:
            bam_rg = add_labs.bam_rg,
            popmapfile = create_popmapFile.popmapfile
    }

    output {
        File refmap = ref_map.stacksVCF
        File freebayes_vcf = freebayes.freebayesVCF
        File gatk_vcf = GenotypeGVCFs.output_vcf
    }
}
