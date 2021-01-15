#!/bin/bash

#SBATCH -J cromwell_main
#SBATCH -c 3
#SBATCH --mem 10240
#SBATCH -p short
#SBATCH -o /data1/aafgarci/cris/main.log
#SBATCH -e /data1/aafgarci/cris/main.err

echo test