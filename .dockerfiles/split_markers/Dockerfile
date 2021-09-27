FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt update \
   && apt-get install -y devscripts software-properties-common

RUN apt update \
    && apt install -y git build-essential vcftools
    
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
    && apt update \
    && apt install -y r-base

RUN R -e "install.packages('ape', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('dplyr', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('magrittr', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('memuse', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('pinfsc50', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('Rcpp', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('stringr', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('tibble', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('vegan', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('viridisLite', repos='http://brieger.esalq.usp.br/CRAN/')"

RUN R -e "packageurl <- 'https://cran.r-project.org/src/contrib/Archive/vcfR/vcfR_1.11.0.tar.gz';\
         install.packages(packageurl, repos=NULL, type='source', dependencies=T)"
         
COPY ./split.R /opt/scripts/
         