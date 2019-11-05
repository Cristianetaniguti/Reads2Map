version 1.0

import "structs/build_pirs_profilesS.wdl"

workflow pirs_profiles{
    input{
        Dataset dataset
    }

    call BaseCalling{
        input:
	        fastq = dataset.fastq_file,
            ref   = dataset.ref,
            ref_idx = dataset.ref_idx
    }

    call GCDepth{
        input:
           file_soap = BaseCalling.file_soap,
           file_single = BaseCalling.file_single,
           ref = dataset.ref
    }

    call Soap2sam {
        input:
            file_soap = BaseCalling.file_soap
    }

    call Indel{
        input:
           sam_file = Soap2sam.sam_file
    }

    output{
        File indel_profile = Indel.indel_profile
        File base_calling_profile = BaseCalling.base_calling_profile
        File GC_bias_100 = GCDepth.GC_bias_100
        File GC_bias_150 = GCDepth.GC_bias_150
        File GC_bias_200 = GCDepth.GC_bias_200
    }
}

task BaseCalling{
    input{
        File fastq
        File ref
        File ref_idx
    }

    command <<<
        soap -a ~{fastq} -D ~{ref_idx} -o align.soap -2 align.single -p 6 -t -s 40 -l 32 -m 600 -x 720 -v 5 -g 0 2> align.log

        ./pirs/baseCalling_Matrix_calculator -r ~{ref} -l 91 -o base_calling -b align.soap
    >>>

    runtime{
	docker:"cristaniguti/soap-pirs"
    }
    
    output{
        File file_soap = "align.soap"
        File file_single = "align.single"
        File base_calling_profile = "base_calling.matrix.gz"
    }

}

task GCDepth{
    input{
        File file_soap
        File file_single
        File ref
    }

    command <<<
        soap.coverage -cvg -onlyuniq -i ~{file_soap} ~{file_single} -refsingle ~{ref} -o align_soap.dresult -depthsingle align_soap.depth > align_soap.deplog 2> align_soap.deperr
        
        ./pirs/gc_coverage_bias -r ~{ref} -o GC_bias -w 100,150,200 align_soap.depth
    >>>

    runtime{
	docker:"cristaniguti/soap-pirs"    
    }

    output{
        File GC_bias_100 = "GC_bias_100.dat"
        File GC_bias_150 = "GC_bias_150.dat"
        File GC_bias_200 = "GC_bias_200.dat"
    }
}


task Soap2sam{
    input{
        File file_soap
    }

    command <<<
        ./BamDeal/bin/BamDeal_Linux convert soap2bam -i ~{file_soap} -s align.sam
    >>>

    runtime{
        docker:"cristaniguti/soap-pirs" 
    }

    output{
        File sam_file = "align.sam"
    }
}

task Indel{
    input{
      File sam_file
    }

    command <<<
        ./pirs/indelstat_sam_bam ~{sam_file} indelstat_profile
    >>>

    runtime{
	    docker:"cristaniguti/soap-pirs"    
    }

    output{
       File indel_profile = "indelstat_profile.InDel.matrix"
    }
}