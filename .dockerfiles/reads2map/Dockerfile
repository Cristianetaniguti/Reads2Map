FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt update \
   && apt-get install -y libcurl4-gnutls-dev \
   libssl-dev \
   cmake \
   libftgl2=2.1.3~rc5-4+nmu1.2 \
   libglu1-mesa-dev=9.0.0-2.1build1 \
   libfreetype6-dev=2.8.1-2ubuntu2.2  \
   gcc-8-base=8.4.0-1ubuntu1~18.04 \
   libgfortran5=8.4.0-1ubuntu1~18.04 \
   python-pip \
   python-tk=2.7.17-1~18.04 \
   libxml2-dev \
   libnlopt-dev=2.4.2+dfsg-4 \
   devscripts=2.17.12ubuntu1.1 \
   software-properties-common=0.96.24.32.22 \
   libgit2-dev=0.26.0+dfsg.1-1.1ubuntu0.2 \
   && pip install matplotlib==2.0.2 numpy 

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' \
    && apt update \
    && apt install -y r-base

RUN R -e 'install.packages("remotes")'
RUN R -e 'remotes::install_github("r-lib/remotes", ref = "97bbf81")'

RUN Rscript -e 'install.packages("R.utils")'
RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'remotes::install_version("stringi",upgrade="never", version = "1.5.3")'
RUN Rscript -e 'remotes::install_version("vcfR",upgrade="never", version = "1.12.0")'
RUN Rscript -e 'remotes::install_version("reshape2",upgrade="never", version = "1.4.4")'
RUN Rscript -e 'remotes::install_version("ggplot2",upgrade="never", version = "3.3.3")'
RUN Rscript -e 'remotes::install_version("doParallel",upgrade="never", version = "1.0.16")'
RUN Rscript -e 'remotes::install_version("matrixStats",upgrade="never", version = "0.58.0")'
RUN Rscript -e 'remotes::install_version("BiocManager",upgrade="never", version = "1.30.12")'
RUN Rscript -e 'remotes::install_version("tidyverse",upgrade="never", version = "1.3.0")'
RUN Rscript -e 'remotes::install_version("largeList",upgrade="never", version = "0.3.0")'
RUN Rscript -e 'remotes::install_version("data.table",upgrade="never", version = "1.14.0")'
RUN Rscript -e 'remotes::install_version("ape",upgrade="never", version = "5.4-1")'
RUN Rscript -e 'remotes::install_version("memuse",upgrade="never", version = "4.1-0")'
RUN Rscript -e 'remotes::install_version("pinfsc50",upgrade="never", version = "1.2.0")'
RUN Rscript -e 'remotes::install_version("vegan",upgrade="never", version = "2.5-7")'
RUN Rscript -e 'remotes::install_version("vroom",upgrade="never", version = "1.4.0")'

RUN Rscript -e 'BiocManager::install("pcaMethods")'
RUN Rscript -e 'BiocManager::install("VariantAnnotation", ask=F)'
RUN Rscript -e 'remotes::install_github("lvclark/polyRAD", ref = "e22ea55372413fe01eb06bb7448f3789a16c8451")'
RUN Rscript -e 'remotes::install_version("gsalib",upgrade="never", version = "2.1")'
RUN Rscript -e 'remotes::install_github("tpbilton/GUSbase", ref = "92119b9c57faa7abeede8236d24a4a8e85fb3df7")'

RUN Rscript -e 'remotes::install_github("tpbilton/GUSMap", ref = "4d7d4057049819d045750d760a45976c8f30dac6")'
RUN Rscript -e 'remotes::install_github("dcgerard/updog")'

RUN Rscript -e 'remotes::install_github("Cristianetaniguti/onemap")'

# Still privates
RUN Rscript -e 'remotes::install_github("Cristianetaniguti/simuscopR")'
RUN Rscript -e 'remotes::install_github("Cristianetaniguti/Reads2MapTools")'
RUN Rscript -e 'remotes::install_github("Cristianetaniguti/MAPpoly" )'
