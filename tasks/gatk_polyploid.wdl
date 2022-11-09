version 1.0

struct Reference {
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File ref_sa
  File ref_amb
  File ref_bwt
  File ref_ann
  File ref_pac
}


workflow GATK_poly {

  input {
    File samples_info
    Reference references
    Int max_cores
    Int chunk_size
    Int ploidy
  }

    call SepareChunks {
        input:
            families_info=samples_info,
            chunk_size = chunk_size
    }

    scatter (chunk in SepareChunks.chunks) {

        Array[Array[String]] sample_file = read_tsv(chunk)

        call RunBwaAlignment {
            input:
                sampleName  = sample_file[1],
                reads       = sample_file[0],
                libraries   = sample_file[2],
                references  = references,
                max_cores   = max_cores,
        }
    }
    
call CreateChunks {
    input:
      bams=flatten(RunBwaAlignment.bam),
      bams_index=flatten(RunBwaAlignment.bai),
      chunk_size=chunk_size,
      reference_fasta = references.ref_fasta
  }

  scatter (chunk in zip(CreateChunks.bams_chunks, CreateChunks.bais_chunks)) {

    call HaplotypeCaller {
      input:
        bams = read_lines(chunk.left),
        bams_index = read_lines(chunk.right),
        reference_fasta = references.ref_fasta,
        reference_fai = references.ref_fasta_index,
        reference_dict = references.ref_dict,
        ploidy = ploidy
    }
  }

  Array[String] calling_intervals = read_lines(CreateChunks.interval_list)

  scatter (interval in calling_intervals) {
    call ImportGVCFs {
      input:
        vcfs=flatten(HaplotypeCaller.vcfs),
        vcfs_index=flatten(HaplotypeCaller.vcfs_index),
        reference_fasta=references.ref_fasta,
        reference_fai=references.ref_fasta_index,
        reference_dict=references.ref_dict,
        interval = interval
    }

    call GenotypeGVCFs {
      input:
        workspace_tar = ImportGVCFs.output_workspace,
        interval = interval,
        reference_fasta = references.ref_fasta,
        reference_fai = references.ref_fasta_index,
        reference_dict = references.ref_dict
    }
  }

  call MergeVCFs {
    input:
      input_vcfs = GenotypeGVCFs.vcf,
      input_vcf_indices = GenotypeGVCFs.vcf_tbi,
      ref_fasta = references.ref_fasta,
      ref_fasta_index = references.ref_fasta_index,
      ref_dict = references.ref_dict
  }

  output {
    File gatk_vcf = MergeVCFs.output_vcf
  }
}


task SepareChunks {
    input {
        File families_info
        Int chunk_size
    }

    command <<<
        R --vanilla --no-save <<RSCRIPT
            df <- read.table("~{families_info}")
            split_df <- split.data.frame(df, df[,2])

            n_chunk <- as.integer(length(split_df)/~{chunk_size})
            chunk_temp <- rep(1:n_chunk, each=~{chunk_size})
            chunk <- c(chunk_temp, rep(n_chunk+1, length(split_df) - length(chunk_temp)))
            chunk_sep <- split(split_df, chunk)

            for(i in 1:length(chunk_sep)){
                df <- do.call(rbind, unlist(chunk_sep[i], recursive = F))
                df <- t(df)
                write.table(df, file = paste0("chunk_",i, ".txt"), quote = F, col.names = F, row.names = F, sep="\t")
            }

        RSCRIPT

    >>>

    runtime {
        job_name: "SepareChunksIndividuals"
        docker: "cristaniguti/r-samtools:latest"
        node:"--nodes=1"
        mem:"--mem=1G"
        tasks:"--ntasks=1"
        time:"00:05:00"
    }

    meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Split the fastq files into chunks to be aligned in parallel in the next task."
    }

    output {
        Array[File] chunks = glob("chunk*")
    }
}


task RunBwaAlignment {

  input {
    Array[String] sampleName
    Array[File] reads
    Array[String] libraries
    Reference references
    Int max_cores
  }

  command <<<
    mkdir tmp
    
    reads_list=( ~{sep=" " reads} )
    lib_list=( ~{sep=" " libraries} )
    sampleName_list=( ~{sep=" " sampleName})
    BAMS=()
    for index in ${!reads_list[*]}; do
      echo "${reads_list[$index]} is in ${lib_list[$index]}"
      bwa_header="@RG\tID:${sampleName_list[$index]}.${lib_list[$index]}\tLB:lib-${lib_list[$index]}\tPL:illumina\tSM:${sampleName_list[$index]}\tPU:FLOWCELL1.LANE1.${lib_list[$index]}"
      /usr/gitc/./bwa mem -t ~{max_cores} -R "${bwa_header}" ~{references.ref_fasta} "${reads_list[$index]}" | \
          java -jar /usr/gitc/picard.jar SortSam \
            I=/dev/stdin \
            O="${sampleName_list[$index]}.${lib_list[$index]}.sorted.bam" \
            TMP_DIR=./tmp \
            SORT_ORDER=coordinate \
            CREATE_INDEX=true;
      mv "${sampleName_list[$index]}.${lib_list[$index]}.sorted.bai" "${sampleName_list[$index]}.${lib_list[$index]}.sorted.bam.bai";
      BAMS+=("I=${sampleName_list[$index]}.${lib_list[$index]}.sorted.bam")
    done

    sampleName_unique=($(echo "${sampleName_list[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

    # Check if there are replicated samples
    for index in ${!sampleName_unique[*]}; do
      NFILES=($(echo ${sampleName_unique[$index]}.*.bam))
      echo ${NFILES[*]}
      echo ${#NFILES[@]}
      REP=()
      if [ "${#NFILES[@]}" -gt 1 ]; then                
        for file in ${!NFILES[*]}; do
          REP+=("I=${NFILES[$file]}")
        done
        echo ${REP[*]}

        java -jar /usr/gitc/picard.jar MergeSamFiles ${REP[*]} \
          O=${sampleName_unique[$index]}.sorted.bam \
          CREATE_INDEX=true \
          TMP_DIR=./tmp
      else
        mv ${sampleName_unique[$index]}.*.bam ${sampleName_unique[$index]}.sorted.bam
        mv ${sampleName_unique[$index]}.*.bai ${sampleName_unique[$index]}.sorted.bam.bai
      fi

    done
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    job_name: "RunBwaAlignment"
    node:"--nodes=1"
    mem:"--mem=32GB"
    tasks:"--ntasks-per-node=11"
    time:"05:00:00"
  }

  meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Run [BWA](http://bio-bwa.sourceforge.net/) MEM alignment."
  }

  output {
    Array[File] bam = glob("*.sorted.bam")
    Array[File] bai = glob("*.sorted.bam.bai")
  }
}


task CreateChunks {
  input {
    Array[String] bams
    Array[String] bams_index
    File reference_fasta
    Int chunk_size
  }

  command <<<
    set -e
    for i in ~{sep=" " bams}; do echo $i >> lof_bams.txt; done
    for i in ~{sep=" " bams_index}; do echo $i >> lof_bais.txt; done

    split -l ~{chunk_size} lof_bams.txt bams.
    split -l ~{chunk_size} lof_bais.txt bais.

    cat ~{reference_fasta} | grep '>' | tr '\n' ',' | sed '$ s/.$//' | sed 's/,/ \n/g' | sed 's/>//g' > intervals.txt  
  >>>

  runtime {
    docker: "ubuntu:20.04"
    job_name: "CreateChunks"
    node:"--nodes=1"
    mem:"--mem=1G"
    cpu:"--ntasks=1"
    time:"00:05:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Split the the samples BAM alignment files in chunks."
  }

  output {
    Array[File] bams_chunks = glob("bams.*")
    Array[File] bais_chunks = glob("bais.*")
    File interval_list = "intervals.txt"
  }
}

task HaplotypeCaller {
  input {
    File reference_fasta
    File reference_dict
    File reference_fai
    Array[File] bams
    Array[File] bams_index
    Int ploidy
  }

  command <<<
    set -euo pipefail

    for bam in ~{sep=" " bams}; do ln -s $bam .; done
    for bai in ~{sep=" " bams_index}; do ln -s $bai .; done

    mkdir vcfs
    ## gvcf for each sample
    for bam in *.bam; do
      out_name=$(basename -s ".bam" "$bam")
      /usr/gitc/gatk4/./gatk --java-options "-Xmx10G -Xms2G" HaplotypeCaller \
        -ERC GVCF \
        -R ~{reference_fasta} \
        -ploidy ~{ploidy} \
        -I "$bam" \
        -O "vcfs/${out_name}.g.vcf.gz" \
        --max-reads-per-alignment-start 0 &
    done

    wait
    
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    job_name: "HaplotypeCaller"
    node:"--nodes=1"
    mem:"--mem=110G" # each sample require 4 cores and 10 GB, here is considering 10 samples/node, 10 extra GB
    tasks:"--ntasks-per-node=42" # 2 extra cores
    time:"24:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Runs [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)."
  }

  output {
    Array[File] vcfs = glob("vcfs/*.vcf.gz")
    Array[File] vcfs_index = glob("vcfs/*.vcf.gz.tbi")
  }
}

task ImportGVCFs  {
  input {
    Array[File] vcfs
    Array[File] vcfs_index
    File reference_fasta
    File reference_fai
    File reference_dict
    String interval
  }

  command <<<
    set -euo pipefail
    grep ">" ~{reference_fasta} | sed 's/^.//' > interval.list
    mkdir gvcfs
    for i in ~{sep=" " vcfs}; do ln -s $i gvcfs/; done
    for i in ~{sep=" " vcfs_index}; do ln -s $i gvcfs/; done

    /usr/gitc/gatk4/./gatk --java-options "-Xmx10G -Xms2G" GenomicsDBImport \
      --batch-size 50 \
      --reader-threads 5 \
      --genomicsdb-workspace-path cohort_db \
      -L ~{interval} \
      -V $(find gvcfs/*.g.vcf.gz -type l | paste -d',' -s | sed 's/,/ -V /g') \
      --consolidate

    tar -cf cohort_db.tar cohort_db

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    job_name: "ImportGVCFs"
    node:"--nodes=1"
    mem:"--mem=20G"
    tasks:"--ntasks-per-node=5"
    time:"24:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Runs [GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport)."
  }

  output {
    File output_workspace = "cohort_db.tar"
  }
}

task GenotypeGVCFs   {
  input {
    File workspace_tar
    File reference_fasta
    File reference_fai
    File reference_dict
    String interval
  }

  command <<<
    set -euo pipefail

    tar -xf ~{workspace_tar}

    /usr/gitc/gatk4/./gatk --java-options "-Xmx10G -Xms2G" GenotypeGVCFs \
      -R ~{reference_fasta} \
      -V gendb://cohort_db \
      -L ~{interval} \
      -G StandardAnnotation \
      -O gatk.vcf.gz 

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    job_name: "GenotypeGVCFs"
    node:"--nodes=1"
    mem:"--mem=20G"
    tasks:"--ntasks=1"
    time:"24:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Runs [GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs)."
  }

  output {
    File vcf = "gatk.vcf.gz"
    File vcf_tbi = "gatk.vcf.gz.tbi"
  }
}

task MergeVCFs {

  input {
    Array[File] input_vcfs
    Array[File] input_vcf_indices
    File ref_fasta
    File ref_fasta_index
    File ref_dict
  }

  command <<<

    /usr/gitc/gatk4/./gatk --java-options "-Xmx10G -Xms2G" \
      MergeVcfs \
        -I ~{sep=' -I' input_vcfs} \
        -O gatk_joint.vcf.gz
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    job_name: "MergeVCFs"
    node:"--nodes=1"
    mem:"--mem=20G"
    tasks:"--ntasks=1"
    time:"10:00:00"
  } 

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Runs [MergeVcfs](https://gatk.broadinstitute.org/hc/en-us/articles/360037226612-MergeVcfs-Picard-)."
  }

  output {
    File output_vcf = "gatk_joint.vcf.gz"
    File output_vcf_index = "gatk_joint.vcf.gz.tbi"
  }
}