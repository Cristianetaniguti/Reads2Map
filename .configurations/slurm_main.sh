#!/bin/bash

#SBATCH -J cromwell_main
#SBATCH -c 3
#SBATCH --mem 10240
#SBATCH -p short
#SBATCH -o /data1/aafgarci/cris/main.log
#SBATCH -e /data1/aafgarci/cris/main.err

java -jar -Dconfig.file=/home/cristiane/github/onemap_workflows/.configurations/cromwell_cache.conf \
     -jar /home/cristiane/cromwell-55.jar \
     run /home/cristiane/github/onemap_workflows/SimulatedReads.wdl \
     -i /home/cristiane/github/onemap_workflows/inputs/SimulatedReads.inputs.toy_sample.json \
     -o /home/cristiane/github/onemap_workflows/.configurations/options.json
