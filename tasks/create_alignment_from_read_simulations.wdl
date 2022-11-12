version 1.0

import "../structs/struct_reads_simu.wdl"
import "custom/alignment.wdl" as alg
import "custom/vcf.wdl"
import "custom/r_libs.wdl"
import "custom/pedigree_simulator_utils.wdl"
import "pirs.wdl"
import "pedigree_simulator.wdl"
import "./utils.wdl" as utils

workflow CreateAlignmentFromSimulation {
    input {
        Reference references
        Family family
        Sequencing sequencing
        Int max_cores
        Int chunk_size
    }

  # User can provide specific variants in VCF file
  # If not, use pirs to simulate
  if (!defined(sequencing.emp_vcf)){
    call pirs.GenerateAlternativeGenome {
      input:
        seed       = family.seed,
        ref_genome = references.ref_fasta
    }

    call pedigree_simulator_utils.CreatePedigreeSimulatorInputs {
      input:
        seed       = family.seed,
        snps       = GenerateAlternativeGenome.snps,
        indels     = GenerateAlternativeGenome.indels,
        cmBymb     = family.cmBymb,
        ref        = references.ref_fasta,
        ref_fai    = references.ref_fasta_index,
        cross      = family.cross,
        popsize    = family.popsize,
        ploidy     = family.ploidy,
        doses      = family.doses
    }
  }

  if (defined(sequencing.emp_vcf)){
    call pedigree_simulator_utils.Vcf2PedigreeSimulator {
      input:
        vcf_file = sequencing.emp_vcf,
        ref_map = sequencing.ref_map,
        seed = family.seed,
        popsize = family.popsize,
        vcf_parent1 = sequencing.vcf_parent1,
        vcf_parent2 = sequencing.vcf_parent2
    }
  }

  File mapfile_sele = select_first([Vcf2PedigreeSimulator.mapfile_map, CreatePedigreeSimulatorInputs.mapfile_nomap])
  File founderfile_sele = select_first([Vcf2PedigreeSimulator.founderfile_map, CreatePedigreeSimulatorInputs.founderfile_nomap])
  File chromfile_sele = select_first([Vcf2PedigreeSimulator.chromfile_map, CreatePedigreeSimulatorInputs.chromfile_nomap])
  File parfile_sele = select_first([Vcf2PedigreeSimulator.parfile_map, CreatePedigreeSimulatorInputs.parfile_nomap])
  File ref_alt_alleles_sele = select_first([Vcf2PedigreeSimulator.ref_alt_alleles_map, CreatePedigreeSimulatorInputs.ref_alt_alleles_nomap])
  File simulated_phases_sele = select_first([Vcf2PedigreeSimulator.simulated_phases_map, CreatePedigreeSimulatorInputs.simulated_phases_nomap])

  call pedigree_simulator.RunPedigreeSimulator {
    input:
      mapfile     = mapfile_sele,
      founderfile = founderfile_sele,
      chromfile   = chromfile_sele,
      parfile     = parfile_sele
  }

  call pedigree_simulator_utils.ConvertPedigreeSimulationToVcf {
    input:
      seed            = family.seed,
      depth           = sequencing.depth,
      genotypes_dat   = RunPedigreeSimulator.genotypes_dat,
      map_file        = mapfile_sele,
      chrom_file      = chromfile_sele,
      ref_alt_alleles = ref_alt_alleles_sele,
      popsize         = family.popsize,
      mapsize         = sequencing.mapsize
  }

  call vcf.GenerateSampleNames {
    input:
      simulated_vcf = ConvertPedigreeSimulationToVcf.simu_vcf
  }

  if(sequencing.library_type == "WGS" || sequencing.library_type == "exome"){

    call r_libs.SimuscopProfile {
      input:
        library_type = sequencing.library_type,
        emp_bam = sequencing.emp_bam,
        vcf     = ConvertPedigreeSimulationToVcf.simu_vcf,
        references = references
    }

    scatter (sampleName in GenerateSampleNames.names) {

      call SimuscopSimulation {
        input:
          library_type  = sequencing.library_type,
          sampleName    = sampleName,
          depth         = sequencing.depth,
          emp_bam       = sequencing.emp_bam,
          vcf           = ConvertPedigreeSimulationToVcf.simu_vcf,
          references    = references,
          chrom         = sequencing.chromosome,
          profile       = SimuscopProfile.profile
      }
    }
  }

    # Two option of RADseq
    # The samples need to be simulated together, otherwise they will be all heterozygous
    if(sequencing.library_type == "sdRAD" || sequencing.library_type == "ddRAD"){
      call RADinitioSimulation{
        input:
          depth          = sequencing.depth,
          depth_parents  = sequencing.depth_parents,
          enzyme1        = sequencing.enzyme1,
          enzyme2        = sequencing.enzyme2,
          simu_vcf       = ConvertPedigreeSimulationToVcf.simu_vcf,
          radinitio_vcf  = ConvertPedigreeSimulationToVcf.radinitio_vcf,
          references     = references,
          pcr_cycles     = sequencing.pcr_cycles,
          insert_size    = sequencing.insert_size,
          insert_size_dev = sequencing.insert_size_dev,
          read_length    = sequencing.read_length,
          library_type   = sequencing.library_type,
          chrom          = sequencing.chromosome,
          names          = GenerateSampleNames.names
      }
    }


  Array[File] fastq = select_first([RADinitioSimulation.fastq_rad, SimuscopSimulation.fastq_seq])

  call SepareChunks {
    input:
      fastqs = fastq,
      chunk_size = chunk_size
  }

  scatter (chunk in SepareChunks.chunks){

    call alg.RunBwaAlignmentSimu {
      input:
        reads      = chunk,
        fastqs     = fastq,
        references = references,
        max_cores  = max_cores,
        rm_dupli   = sequencing.rm_dupli 
    }
  }

    # Store for MCHap 
    call utils.MergeBams {
        input:
            bam_files = flatten(RunBwaAlignmentSimu.bam)
    }

  output {
      Array[File] bam = flatten(RunBwaAlignmentSimu.bam)
      Array[File] bai = flatten(RunBwaAlignmentSimu.bai)
      File ref_alt_alleles = ref_alt_alleles_sele
      Array[String] names = GenerateSampleNames.names
      File true_vcf = ConvertPedigreeSimulationToVcf.simu_vcf
      File simu_haplo = ConvertPedigreeSimulationToVcf.simu_haplo
      File simulated_phases = simulated_phases_sele
      File merged_bam = MergeBams.merged_bam
  }

}


task SimuscopSimulation{
 input {
    String library_type
    String sampleName
    Int depth
    File? emp_bam
    File vcf
    Reference references
    String chrom
    File profile
  }

  Int disk_size = ceil(size(emp_bam, "GiB") + size(vcf, "GiB") + size(references.ref_fasta, "GiB") * depth) 
  Int memory_size = 10000

  command <<<
    R --vanilla --no-save <<RSCRIPT
      vcfR.object <- read.vcfR("~{vcf}")

      variants <- vcf2variants(vcfR.object, sample = "~{sampleName}", chrom = "~{chrom}")

      write.table(variants$SNVs, file = "SNVs.txt", sep = "\t", quote = F, col.names = F, row.names = F)
      write.table(variants$indels, file = "indels.txt", sep = "\t", quote = F, col.names = F, row.names = F)
      write.table(variants$insertions, file = "insertions.txt", sep = "\t", quote = F, col.names = F, row.names = F)

      system("cat SNVs.txt indels.txt insertions.txt > variants.txt")

      if("~{library_type}" == "exome"){
        simuReads(ref = "~{references.ref_fasta}",
              profile = "~{profile}",
              variation = "variants.txt",
              target = "bed_file",
              name = "~{sampleName}",
              output = ".",
              layout = "SE",
              threads = 6,
              verbose = 1,
              coverage = ~{depth})
      } else {
        simuReads(ref = "~{references.ref_fasta}",
              profile = "profile",
              variation = "variants.txt",
              name = "~{sampleName}",
              output = ".",
              layout = "SE", # only single-end by now
              threads = 6,
              verbose = 1,
              coverage = ~{depth})
      }

    RSCRIPT

  >>>

  runtime {
    docker: "cristaniguti/reads2map:0.0.1"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "SimuscopSimulation"
    mem:"~{memory_size}M"
    time:"10:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Run [simuReads](https://github.com/qasimyu/simuscop) to simulated exome or WGS sequencing reads."
  }

  output {
    File fastq_seq = "~{sampleName}.fq"
  }
}


task RADinitioSimulation{
  input {
    File simu_vcf
    File radinitio_vcf
    String enzyme1
    String? enzyme2
    Reference references
    Int depth
    Int depth_parents
    Int? insert_size
    Int? insert_size_dev
    Int? pcr_cycles
    Int? read_length
    String library_type
    Array[String] names
    String chrom
  }

  # Difficult to guess how much disk we'll need here.
  Int disk_size = ceil(size(simu_vcf, "GB") + size(radinitio_vcf, "GB") + size(references.ref_fasta, "GB") + 5)
  Int memory_size = 10000  

  command <<<

    echo -e ~{sep=" " names} > temp
    echo "~{chrom}" > chrom.list

    tr -s ' ' '\n' < temp > temp2
    sed 's/$/\tpop0/' temp2 > popmap.tsv

    mkdir simu_inputs_progeny simu_inputs_parents \
          results_progeny results_parents

    mkdir simu_inputs_progeny/msprime_vcfs simu_inputs_progeny/ref_loci_vars \
          simu_inputs_parents/msprime_vcfs simu_inputs_parents/ref_loci_vars

    # Separate progeny from parents because of the different depths
    head -n 2 popmap.tsv > simu_inputs_parents/popmap.tsv
    lines=$(wc -l popmap.tsv | cut -f1 -d' ')
    tail -n $((lines -2)) popmap.tsv > simu_inputs_progeny/popmap.tsv

    vcftools --vcf ~{simu_vcf} --indv P1 --indv P2 --recode --out parents
    vcftools --vcf ~{simu_vcf} --remove-indv P1 --remove-indv P2 --recode --out progeny

    vcftools --vcf ~{radinitio_vcf} --indv P1 --indv P2 --recode --out parents.rad
    vcftools --vcf ~{radinitio_vcf} --remove-indv P1 --remove-indv P2 --recode --out progeny.rad

    gzip parents.recode.vcf
    gzip parents.rad.recode.vcf
    mv parents.rad.recode.vcf.gz simu_inputs_parents/msprime_vcfs/~{chrom}.vcf.gz
    mv parents.recode.vcf.gz simu_inputs_parents/ref_loci_vars/ri_master.vcf.gz

    gzip progeny.recode.vcf
    gzip progeny.rad.recode.vcf
    mv progeny.rad.recode.vcf.gz simu_inputs_progeny/msprime_vcfs/~{chrom}.vcf.gz
    mv progeny.recode.vcf.gz simu_inputs_progeny/ref_loci_vars/ri_master.vcf.gz


    # progeny
    radinitio --make-library-seq \
              --genome ~{references.ref_fasta} \
              --chromosomes chrom.list \
              --out-dir results_progeny/ \
              --make-pop-sim-dir simu_inputs_progeny/ \
              --library-type ~{library_type} \
              --enz ~{enzyme1} \
              ~{"--enz2 " + enzyme2} \
              --insert-mean ~{default="350" insert_size} \
              --insert-stdev ~{default="35" insert_size_dev} \
              --pcr-cycles ~{default="9" pcr_cycles} \
              --coverage ~{default="20" depth} \
              --read-length ~{default="150" read_length}

    # parents
    radinitio --make-library-seq \
          --genome ~{references.ref_fasta} \
          --chromosomes chrom.list \
          --out-dir results_parents/ \
          --make-pop-sim-dir simu_inputs_parents/ \
          --library-type ~{library_type} \
          --enz ~{enzyme1} \
          ~{"--enz2 " + enzyme2} \
          --insert-mean ~{default="350" insert_size} \
          --insert-stdev ~{default="35" insert_size_dev} \
          --pcr-cycles ~{default="9" pcr_cycles} \
          --coverage ~{default="20" depth_parents} \
          --read-length ~{default="150" read_length}

    # Add fake phred score of 40 (H in Illumina 1.8+ Phred+33)
    # Only in forward read
    for i in results_progeny/rad_reads/*.1.fa.gz; do /seqtk/./seqtk seq -F 'I' $i > $(basename ${i/.fa.gz}.fq); done
    for i in results_parents/rad_reads/*.1.fa.gz; do /seqtk/./seqtk  seq -F 'I' $i > $(basename ${i/.fa.gz}.fq); done

  >>>

  runtime{
    docker: "cristaniguti/radinitio:0.0.1"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    maxRetries: 5
    preemptible: 3
    # Slurm
    job_name: "RADinitioSimulation"
    mem:"~{memory_size}M"
    time:"05:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Run [RARinitio](https://github.com/qasimyu/simuscop) to simulated RADseq sequencing reads."
  }

  output {
    Array[File] fastq_rad = glob("*.fq")
  }
}


task SepareChunks {
  input {
      Array[File] fastqs
      Int chunk_size
  }

  Int disk_size = ceil(size(fastqs, "GiB") * 2) 
  Int memory_size = 1000

  command <<<
        R --vanilla --no-save <<RSCRIPT

            files <- c("~{sep="," fastqs}")
            files <- unlist(strsplit(files, split = ","))

            n_chunk <- floor(length(files)/~{chunk_size})

            chunk_temp <- rep(1:n_chunk, each = ~{chunk_size})
            chunk <- c(chunk_temp, rep(n_chunk+1, length(files) - length(chunk_temp)))

            chunk_sep <- split(files, chunk)

            for(i in 1:length(chunk_sep)){
              write.table(chunk_sep[[i]], file = paste0("chunk_",i, ".txt"), quote = F, col.names = F, row.names = F, sep="\t")
            }

        RSCRIPT

  >>>

  runtime {
      docker: "cristaniguti/reads2map:0.0.1"
      cpu:1
      # Cloud
      memory:"~{memory_size} MiB"
      disks:"local-disk " + disk_size + " HDD"
      # Slurm
      job_name: "SepareChunksIndividuals"
      mem:"~{memory_size}M"
      time:"00:10:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Split the simulated fastq files into chunks to be aligned in parallel in the next task."
  }

  output {
    Array[File] chunks = glob("chunk*")
  }
}