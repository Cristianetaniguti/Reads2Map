#!/bin/bash

#SBATCH --export=NONE
#SBATCH --job-name=cromwell_main2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=5G
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chtaniguti@tamu.edu
#SBATCH -o /scratch/user/chtaniguti/cromwell_main2.log
#SBATCH -e /scratch/user/chtaniguti/cromwell_main2.err

ml Java/13.0.2
ml WebProxy

export SINGULARITY_CACHEDIR=/scratch/user/chtaniguti/.singularity
export SINGULARITY_TMPDIR=/scratch/user/chtaniguti/temp

for i in /scratch/user/chtaniguti/euc_raw/A_1/A_1_raw/*.fastq.gz; do # Will run in all files in the directory finishing with sub.fastq
    echo $i
    singularity run --bind /scratch/user/chtaniguti/euc_raw/A_1/A_1_raw/:/opt \
     /scratch/user/chtaniguti/.singularity/biocontainers_fastqc_v0.11.9_cv7.sif fastqc /opt/$i
done

singularity run --bind /scratch/user/chtaniguti/euc_raw/A_1/A_1_raw/:/scratch/user/chtaniguti/euc_raw/A_1/A_1_raw/ \
 -W /scratch/user/chtaniguti/euc_raw/A_1/A_1_raw/ /scratch/user/chtaniguti/.singularity/ewels_multiqc.sif . --title "A_1"


# TODO: will we keep it?
