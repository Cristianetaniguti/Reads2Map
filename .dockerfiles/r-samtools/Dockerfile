FROM r-base:3.6.0

RUN apt update \
   && apt-get install -y samtools


COPY ./vcf2pedigreeSim.R /opt/scripts/