task createaltgenome{
     File ref_genome

     command{
	pirs diploid ${ref_genome} -s 0.001 -d 0 -v 0 -o alt
     }
     runtime{
	docker:"wkusmirek/pirs:latest"
	memory: "1 GB"
        cpu: "1"
     }
     output{
	File alt_fasta = "alt.snp.fa"
     	File snps = "alt.snp.lst"
     }
}

workflow creategenome{
	 File ref
	 
	 call createaltgenome{
	      input:ref_genome=ref 
	 }
}