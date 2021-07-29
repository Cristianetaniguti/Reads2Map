#!/bin/bash

#SBATCH --export=NONE
#SBATCH -J cromwell_main
#SBATCH --nodes=1                    
#SBATCH --mem=1G
#SBATCH --time=01:30:00
#SBATCH -o /scratch/user/chtaniguti/main.log
#SBATCH -e /scratch/user/chtaniguti/main.err

java -jar -Dconfig.file=/scratch/user/chtaniguti/Reads2Map/.configurations/cromwell_cache.conf \
     -jar /scratch/user/chtaniguti/cromwell-65.jar \
     run /scratch/user/chtaniguti/Reads2Map/tasks/snpcalling_emp.wdl \
     -i /scratch/user/chtaniguti/Reads2Map/inputs/toy_sample_emp.inputs.json 
