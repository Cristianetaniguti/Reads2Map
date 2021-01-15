#!/bin/bash

#SBATCH -J cromwell_main
#SBATCH -c 3
#SBATCH --mem 10240
#SBATCH -p short
#SBATCH -o /data1/aafgarci/cris/main.log
#SBATCH -e /data1/aafgarci/cris/main.err

java -jar -Dconfig.file=/data1/aafgarci/cris/simulations/onemap_workflows/.configurations/cromwell_sing.conf \
     -jar /data1/aafgarci/cris/cromwell-55.jar \
     run /data1/aafgarci/cris/simulations/onemap_workflows/SimulatedReads.wdl \
     -i /data1/aafgarci/cris/simulations/onemap_workflows/inputs/SimulatedReads.depth20.popsize200.multi.inputs 