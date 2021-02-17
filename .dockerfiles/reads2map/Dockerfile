FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt update \
   && apt-get install -y libcurl4-gnutls-dev \ 
   libssl-dev=1.1.1-1ubuntu2.1~18.04.7 \
   libftgl2=2.1.3~rc5-4+nmu1.2 \
   libglu1-mesa-dev=9.0.0-2.1build1 \
   libfreetype6-dev=2.8.1-2ubuntu2.1  \
   gcc-8-base=8.4.0-1ubuntu1~18.04 \
   libgfortran5=8.4.0-1ubuntu1~18.04 \
   python-pip=9.0.1-2.3~ubuntu1.18.04.4 \
   python-tk=2.7.17-1~18.04 \
   libxml2-dev=2.9.4+dfsg1-6.1ubuntu1.3 \
   libnlopt-dev=2.4.2+dfsg-4 \
   devscripts=2.17.12ubuntu1.1 \
   software-properties-common=0.96.24.32.14 \
   libgit2-dev=0.26.0+dfsg.1-1.1ubuntu0.2 \
   && pip install matplotlib==2.0.2 numpy

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' \
    && apt update \
    && apt install -y r-base

RUN R -e "install.packages('R.utils', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('devtools', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('stringi', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('vcfR', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('reshape2', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('ggplot2', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('doParallel', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('matrixStats', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('BiocManager', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('tidyverse', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('largeList', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('data.table', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('ape', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('memuse', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('pinfsc50', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('vegan', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('vroom', repos='http://brieger.esalq.usp.br/CRAN/');\
          BiocManager::install('pcaMethods');\
          BiocManager::install('VariantAnnotation', ask=F);\
          install.packages('polyRAD', repos='http://brieger.esalq.usp.br/CRAN/');\
          devtools::install_github('tpbilton/GUSbase');\
          devtools::install_github('tpbilton/GUSMap')"

RUN R -e "devtools::install_github('Cristianetaniguti/onemap' )"

# Still in private repo
COPY ./genotyping4onemap_0.1.7.tar.gz /opt/
RUN R -e "install.packages('/opt/genotyping4onemap_0.1.7.tar.gz', type = 'source', repos=NULL)"

COPY ./functions_simu.R /opt/scripts/

COPY ./functions_empirical.R /opt/scripts
