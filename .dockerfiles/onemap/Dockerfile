FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y \
        libcurl4-gnutls-dev \
        libssl-dev \
        libftgl2 \
        libglu1-mesa-dev \
        libftgl2 \
        libfreetype6-dev  \
        gcc-8-base \
        libgfortran5 \
        python-pip \
        python-tk \
        libxml2-dev \
        libnlopt-dev \
        devscripts \
        software-properties-common \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
    && apt-get update \
    && apt-get install -y \
        r-base \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

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
          BiocManager::install('pcaMethods');\
          install.packages('polyRAD', repos='http://brieger.esalq.usp.br/CRAN/');\
          devtools::install_github('tpbilton/GUSbase');\
          devtools::install_github('tpbilton/GUSMap@v2.0.0')"

RUN R -e "packageurl <- 'https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.3.tar.gz';\
          install.packages(packageurl, repos=NULL, type='source');\
          packageurl <- 'https://cran.r-project.org/src/contrib/Archive/vcfR/vcfR_1.11.0.tar.gz';\
          install.packages(packageurl, repos=NULL, type='source');\
          BiocManager::install('VariantAnnotation', ask=F)"

RUN apt-get update \
    && pip install matplotlib==2.0.2 numpy

RUN R -e "devtools::install_github('Cristianetaniguti/onemap')"

COPY ./functions_simu.R /opt/scripts/

COPY ./functions_empirical.R /opt/scripts

COPY ./genotyping4onemap_0.1.7.tar.gz /opt/
RUN R -e "install.packages('/opt/genotyping4onemap_0.1.7.tar.gz', type = 'source', repos=NULL)"

RUN R -e "packageurl <- 'https://cran.r-project.org/src/contrib/Archive/updog/updog_1.2.0.tar.gz';\
          install.packages(packageurl, repos=NULL, type='source')"

