task create_alt_genome{
	File ref_genome

     command{
		/pirs/src/pirs/./pirs diploid ${ref_genome} -s 0.001 -d 0 -v 0 -o alt
		chmod 777 "alt.snp.lst" "alt.snp.fa"
     }
     runtime{
		docker:"pirs"
     }
     output{
		File alt_fasta = "alt.snp.fa"
    		File snps = "alt.snp.lst"
     }
}

task pedsim_files{
	File snp_file
	String genome_size
	String cmBymb
	File R_script
	command{
		Rscript --vanilla ${R_script} ${snp_file} ${genome_size}  ${cmBymb}
		chmod 777 "mapfile.map" "founderfile.gen" "sim.par" "inb.chrom"
	}
	runtime{
		docker:"r-base:3.6.0"
	}
	output{
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
	File par_file
	File chrom_file
	command{
		java -jar ${pedigreeSimJar} ${par_file} 
	}
	output{
		File genotypes_dat = "sim_inb_genotypes.dat"
	}
}

workflow F2{
	 File ref
	 String genome_size
	 String cmBymb
	 File R_script
	 File pedigreeSim_jar
	 
	 call create_alt_genome{
	     input:ref_genome=ref 
	 }
	 call pedsim_files{
		 input:snp_file = create_alt_genome.snps,
		 genome_size=genome_size,
		 cmBymb=cmBymb,
		 R_script=R_script
	 }
	  call pedigreeSim{
	 	 input: map_file = pedsim_files.mapfile,
	 	 founder_file= pedsim_files.founderfile, 
	 	 par_file=pedsim_files.parfile,
	 	 chrom_file=pedsim_files.chromfile,
	 	 pedigreeSimJar=pedigreeSim_jar
	 }
}