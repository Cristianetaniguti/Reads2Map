FROM r-base:4.0.0

RUN apt update \
   && apt-get install -y samtools

RUN R -e "install.packages('vcfR', repos='http://brieger.esalq.usp.br/CRAN/')"

COPY ./onemapUTILS_0.2.0.tar.gz /opt/
RUN R -e "install.packages('/opt/onemapUTILS_0.2.0.tar.gz', type = 'source', repos=NULL)"