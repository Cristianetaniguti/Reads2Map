task create_alt_genome{
	File ref_genome

     command{
		/pirs/src/pirs/./pirs diploid ${ref_genome} -s 0.001 -d 0 -v 0 -o alt
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
	File mapfile
	File founderfile
	File parfile
	File chromfile
	
	command{
		PATH=$PATH:${mapfile}
		PATH=$PATH:${chromfile}
		PATH=$PATH:${founderfile}
		java -jar ${pedigreeSimJar} ${parfile}
	}
	runtime{
		docker:"java"
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
	 	 input: mapfile = pedsim_files.mapfile,
	 	 founderfile= pedsim_files.founderfile, 
	 	 parfile=pedsim_files.parfile,
	 	 chromfile=pedsim_files.chromfile,
	 	 pedigreeSimJar=pedigreeSim_jar
	 }
}