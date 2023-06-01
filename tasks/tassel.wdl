version 1.0


task transposeSamples {
  input {
    File families_info
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT

      file <- read.table("~{families_info}")
      file <- t(file)
      write.table(file, file = "transpose_sample.tsv", quote = F, row.names = F, col.names = F, sep = "\t")

    RSCRIPT
  >>> 

  runtime {
    docker: "cristaniguti/r-samtools:latest"
    singularity: "docker://cristaniguti/r-samtools:latest"
    cpu: 1
    # Cloud
    memory:"1 GiB"
    disks:"local-disk 3 GiB HDD"
    # Slurm
    job_name: "transposeSamples"
    mem:"1G"
    time: 4
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "transpose the sample info"
  }

  output {
    File transpose_samples = "transpose_sample.tsv"
  }  
}

task BarcodeFaker {
  input {
    Array[File] fastq
    Array[String] FullSampleName
  }

  command <<<
    R --vanilla --no-save <<RSCRIPT
      # Marlee's functions
      #this function recursively makes a vector of unique fake barcodes to ensure no non-unique barcodes are generated
      #this step does not protect against enzyme cut site being found in the barcode
      barcode_inventor <- function(fake_barcodes, r, fastqs){
        fake_barcodes <- c()
        
        # make some barcodes of the length needed to ensure unique barcodes for all files are possible
        for(i in 1:length(fastqs)){
          fake_barcodes[i] <- paste(sample(c("A", "T", "C", "G"), size = r, replace = TRUE), collapse = "")
        }
        
        # check if the barcodes are all unique
        if(length(fake_barcodes) == length(unique(fake_barcodes))){
          return(fake_barcodes)
        } else {
          barcode_inventor
        }
      }
      
      # This function adds fake barcodes to the start of reads in Illumina FASTQ files for compatibility with TASSEL5 GBSv2.
      # fastq_dir is the directory path (string) of *unzipped* FASTQ files
      # Returns unzipped FASTQ files with fake barcodes and quality scores, named according to the GBSv2 requirement,
      # and a key, "barcode_key.csv", which matches original file name, new file name, flowcell, lane, and barcode
      # All files are output to the working directory.
      barcode_faker <- function(fastq_dir, FullSampleName){
          print("Always run barcode_faker on a folder containing all files you intend to use in a given GBSv2 run. Multiple separate runs of barcode_faker can cause barcodes to be repeated across sample files, though unlikely, and will certainly cause file names to be repeated.")
          
          
          #get the fastq file paths
          fastqs <- list.files(path = fastq_dir, full.names = TRUE)
          fastqs <- fastqs[which(grepl(".fastq", fastqs) | grepl(".fq", fastqs))]
          fastq.check <- list.files(path = fastq_dir, full.names = FALSE)
          fastq.check <- fastq.check[which(grepl(".fastq", fastq.check) | grepl(".fq", fastq.check))]
          if(length(fastqs) == 0) stop("Files not fount.")
          
          #check that all files in dir are fastqs         #bug fixed 28 July 2020- reported by Dr. Smit Dhakal
          long <- length(grep(".fastq", fastq.check))
          short <- length(grep(".fq", fastq.check))
          if(long != 0){
            if(length(fastq.check) != long){
              stop("Not all of the files in your input folder seem to have .fq or .fastq extensions. Please move or delete such files.")
            }
          }
          if(short != 0){
            if(length(fastq.check) != short){
              stop("Not all of the files in your input folder seem to have .fq or .fastq extensions. Please move or delete such files.")
            }
          }
          
          
          # new method: given a number of input files, determine the min barcode length needed to ensure enough unique barcodes are generated
          # the old way failed if a lot of files had short barcodes of the same length- reported 18 Aug 2022 by Meghan Brady
          # https://www.calculatorsoup.com/calculators/discretemathematics/combinationsreplacement.php
          # n = 4
          # r = min barcode length
          # C = number of files
          # solve for r considering n and C are given
          C <- length(fastqs)
          constant <- (factorial(3) * C) - 6 # get the constant
          r <- round(polyroot(c(constant, -11, -6, -1)), 2) # get the roots of the polynomial
          r <- ceiling(Re(r)[Re(r) > 0]) # subset the real and positive roots, and round up
          if(r < 4){r <- 4} # keep a minimum barcode length of 4, probably unnecessary

          fake_barcodes <- barcode_inventor(fake_barcodes, r, fastqs)
          
          #make the fake perfect quality scores matching the barcode length
          fake_quals <- rep(paste(rep("E", times = r), collapse = ""), times = length(fastqs))
                    
          
          #make fake header files for FASTQ reads
          #headers are totally fake (do not refer to input file)
          #fake headers allows interoperability between TASSEL and demultiplexing software FASTQ formats
          fake_flowcell_no <- seq(from = 1, to = length(fastqs), by = 1)
          fake_headers <- paste("D00553R:56:", "C8B56ANXX", fake_flowcell_no, ":3:1101:1203:2037 1:N:0:3", sep = "")
          fake_flowcells <- paste("C8B56ANXX", fake_flowcell_no, sep = "")
          fake_lanes <- rep(3, times = length(fastqs))
          
          #make fake file names that appropriately reference fake flowcell and lane
          newfile_names <- paste(fake_flowcells, "_", "3","_","fastq.fastq", sep = "")
          
          #stop if files with the fake file names already exist
          workdir <- getwd()
          dirfil <- list.files(path = workdir, full.names = FALSE)
          if(sum(newfile_names %in% dirfil) != 0){
            stop("Files in your working directory seem to have same names as those
                output by this function. Please move them, delete them, or choose a new working directory.")
          }
          
          #write a key of barcode, flowcell, lane and new file names
          key <- cbind(Flowcell = fake_flowcells, Lane = fake_lanes, Barcode = fake_barcodes, FullSampleName)
          write.table(key, "fakebarcodes_key.txt", row.names = FALSE, sep = "\t", quote = FALSE)
          
          #read 400 lines per file at a time, paste on the barcodes and quality scores, make unique flowcells, 
          #and write them out with unique names
          for(i in 1:length(fastqs)){
            incon <- file(fastqs[i], open = "r")
            outcon <- file(newfile_names[i], open = "w")
            while(length(mylines <- readLines(incon, n = 400, warn = FALSE))){
              head_pos <- seq(1, length(mylines), by = 4)
              seq_pos <- seq(2, length(mylines), by = 4)
              qual_pos <- seq(4, length(mylines), by = 4)
              
              for(j in 1:length(mylines)){
                mylines[head_pos[j]] <- fake_headers[i]
                mylines[seq_pos[j]] <- paste(fake_barcodes[i], mylines[seq_pos[j]], sep = "")
                mylines[qual_pos[j]] <- paste(fake_quals[i], mylines[qual_pos[j]], sep = "")
              }
              writeLines(mylines, outcon)
            }
            print(paste("Finished writing File", i, "...", sep = " "))
            #close the connections
            close(outcon, warn = FALSE)
            close(incon, warn = FALSE)
          }
      }

     file_names <- "~{sep=',' fastq}"
     file_names <- unlist(strsplit(file_names, ","))
     is_gz <- basename(file_names[1])
     if(grepl("gz", is_gz)) system(paste("gunzip", file_names))
     dir_name <- dirname(file_names[1])
     
     sample_names <- "~{sep=',' FullSampleName}"
     sample_names <- unlist(strsplit(sample_names, ","))

     barcode_faker(fastq_dir = dir_name, FullSampleName = sample_names)

    RSCRIPT
  >>>

  runtime {
    docker: "cristaniguti/r-samtools:latest"
    singularity: "docker://cristaniguti/r-samtools:latest"
    cpu: 1
    # Cloud
    memory:"2 GiB"
    disks:"local-disk 3 GiB HDD"
    # Slurm
    job_name: "BarcodeFaker"
    mem:"2G"
    time: 4
  }

  meta {
    author: "Marlee Labroo function adapted by Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Add fake barcode adaptors to fastq sequences to be input for TASSEL"
  }

  output {
    File key_file = "fakebarcodes_key.txt"
    Array[File] barcode_fastq = glob("*.fastq")
  }
}

task TasselBeforeAlign {
    input {
        Array[File] fastq
        String enzyme = "ApeKI"
        File key_file
        Int max_ram = 10000        
    }

    Int disk_size = ceil(size(fastq, "GiB"))
    Int memory_min = ceil(max_ram/2)

  command <<<

    mkdir fastqs
    ln -s ~{sep=" " fastq} fastqs

    /usr/tassel/run_pipeline.pl -Xms~{memory_min}m -Xmx~{max_ram}m -fork1 -GBSSeqToTagDBPlugin -e ~{enzyme} -i fastqs \
                             -db GBSV2.db \
                             -k ~{key_file} \
                             -kmerLength 64 -minKmerL 20 -mnQS 20 -mxKmerNum 100000000 -endPlugin -runfork1

    /usr/tassel/run_pipeline.pl -Xms~{memory_min}m -Xmx~{max_ram}m -fork1 -TagExportToFastqPlugin -db GBSV2.db \
      -o tagsForAlign.fa.gz -c 1 -endPlugin -runfork1

  >>>

  runtime {
    docker: "cristaniguti/java-in-the-cloud:0.0.2"
    singularity: "docker://cristaniguti/java-in-the-cloud:0.0.2"
    cpu: 1
    # Cloud
    memory:"~{max_ram} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "GBSSeqToTagDBPlugin"
    mem:"~{max_ram}M"
    time: 24
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "GBSSeqToTagDBPlugin takes fastQ files as input, identifies tags and the taxa in which they appear, and stores this data to a local database. Retrieves distinct tags stored in the database and reformats them to a FASTQ file"
  }

  output {
    File tassel_database = "GBSV2.db"
    File fastq_align = "tagsForAlign.fa.gz"
  }
}

task TasselAfterAlign {
    input {
      File tassel_database
      File bam
      Int max_ram = 10000
      String enzyme
      File key_file
      Array[File] fastq        
    }

    Int disk_size = ceil(size(tassel_database, "GiB"))
    Int memory_min = ceil(max_ram/2)

    command <<<

      samtools view -h ~{bam} > file.sam 
      mv ~{tassel_database} .

      /usr/tassel/run_pipeline.pl -Xms~{memory_min}m -Xmx~{max_ram}m -fork1 -SAMToGBSdbPlugin \
      -i file.sam  \
      -db GBSV2.db \
      -aProp 0.0 -aLen 0 -endPlugin -runfork1

      /usr/tassel/run_pipeline.pl -Xms~{memory_min}m -Xmx~{max_ram}m -fork1 -DiscoverySNPCallerPluginV2 \
      -db GBSV2.db \
      -mnLCov 0.1 -mnMAF 0.01 -deleteOldData true -endPlugin -runfork1

      /usr/tassel/run_pipeline.pl -Xms~{memory_min}m -Xmx~{max_ram}m  -fork1 -SNPQualityProfilerPlugin \
      -db GBSV2.db \
      -deleteOldData true -endPlugin -runfork1

      mkdir fastqs
      ln -s ~{sep=" " fastq} fastqs

      /usr/tassel/run_pipeline.pl -Xms~{memory_min}m -Xmx~{max_ram}m -fork1 -ProductionSNPCallerPluginV2 \
      -db GBSV2.db \
      -e ~{enzyme} -i fastqs \
      -k ~{key_file} \
      -kmerLength 64 \
      -o tassel.vcf -endPlugin -runfork1

    >>>

    runtime {
      docker: "cristaniguti/java-in-the-cloud:0.0.2"
      singularity: "docker://cristaniguti/java-in-the-cloud:0.0.2"
      cpu: 1
      # Cloud
      memory:"~{max_ram} MiB"
      disks:"local-disk " + disk_size + " HDD"
      # Slurm
      job_name: "SAMToGBSdbPlugin"
      mem:"~{max_ram}M"
      time: 24
    }

    meta {
      author: "Cristiane Taniguti"
      email: "chtaniguti@tamu.edu"
      description: "Reads a SAM file to determine the potential positions of Tags against the reference genome. Identifies SNPs from the aligned tags"
    }

    output {
      File tassel_vcf = "tassel.vcf"
    }

}
