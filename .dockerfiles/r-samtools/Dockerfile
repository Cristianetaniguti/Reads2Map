FROM r-base:4.0.0

RUN apt update \
   && apt-get install -y samtools

RUN R -e "install.packages('vcfR', repos='http://brieger.esalq.usp.br/CRAN/')"

COPY ./vcf2pedigreeSim.R /opt/scripts/