version 1.0

import "../structs/dna_seq_structs.wdl"

task RADinitioSimulation {
  input {
    File simu_vcf
    File radinitio_vcf
    String enzyme1
    String? enzyme2
    ReferenceFasta references
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

  runtime {
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
