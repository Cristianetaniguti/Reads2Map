version 1.0

import "structs/build_pirs_profilesS.wdl"

workflow pirs_profiles{
    input{
        Dataset dataset
    }

    call SoapAlign{
        input:
            fastq = dataset.fastq_file,
            ref   = dataset.ref,
            ref_amb = dataset.ref_amb,
            ref_fmv = dataset.ref_fmv,
            ref_pac = dataset.ref_pac,
            ref_sa = dataset.ref_sa,
            ref_sai = dataset.ref_sai,
            ref_ann = dataset.ref_ann,
            ref_hot = dataset.ref_hot,
            ref_bwt = dataset.ref_bwt,
            ref_lkt = dataset.ref_lkt,
            ref_rev_pac = dataset.ref_rev_pac,
            ref_rev_bwt = dataset.ref_rev_bwt,
            ref_rev_fmv = dataset.ref_rev_fmv,
            ref_rev_lkt = dataset.ref_rev_lkt,
            ref_rev_pac = dataset.ref_rev_pac
    }

    call GCDepth{
        input:
           file_soap = SoapAlign.file_soap,
           ref = dataset.ref
    }

    call Soap2sam {
        input:
            file_soap = SoapAlign.file_soap
    }

    call BaseCalling{
        input:
            ref = dataset.ref,
            sam_file = Soap2sam.sam_file
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


task SoapAlign{
    input{
        File fastq
        File ref
        File ref_amb 
        File ref_fmv
        File ref_pac
        File ref_sa
        File ref_sai
        File ref_ann 
        File ref_hot
        File ref_bwt
        File ref_lkt
        File ref_rev_pac
        File ref_rev_bwt
        File ref_rev_fmv
        File ref_rev_lkt
        File ref_rev_pac
    }

    command <<<
        /opt/conda/bin/soap -a ~{fastq} -D ~{ref}.index -o align.soap  -p 10 -t -s 40 -l 32  -v 5 -g 0 2> align.log
    >>>

    runtime{
	docker:"cristaniguti/soap-pirs"
    }
    
    output{
        File file_soap = "align.soap"
    }

}

task GCDepth{
    input{
        File file_soap
        File ref
    }

    command <<<
         /opt/conda/bin/soap.coverage -cvg -onlyuniq -i ~{file_soap}  -refsingle ~{ref} -o align_soap.dresult -depthsingle align_soap.depth > align_soap.deplog 2> align_soap.deperr
        
        /pirs/gc_coverage_bias -r ~{ref} -o GC_bias -w 100,150,200 align_soap.depth
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
        /BamDeal/bin/BamDeal_Linux convert soap2bam -i ~{file_soap} -s align.sam
    >>>

    runtime{
        docker:"cristaniguti/soap-pirs" 
    }

    output{
        File sam_file = "align.sam.gz"
    }
}

task Indel{
    input{
      File sam_file
    }

    command <<<
        /pirs/indelstat_sam_bam ~{sam_file} indelstat_profile
    >>>

    runtime{
	    docker:"cristaniguti/soap-pirs"    
    }

    output{
       File indel_profile = "indelstat_profile.InDel.matrix"
    }
}


task BaseCalling {
    input{
      File sam_file
      File ref
    }

    command <<<
      /pirs/baseCalling_Matrix_calculator -r ~{ref} -l 91 -o base_calling -i ~{sam_file}
    >>>

   runtime{
      docker:"cristaniguti/soap-pirs"
   }

   output{
     File base_calling_profile = "base_calling.matrix.gz"
   }
}
