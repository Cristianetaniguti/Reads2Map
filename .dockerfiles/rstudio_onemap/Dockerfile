FROM cristaniguti/onemap_git:latest

# For supermassa
RUN apt-get update \
    && apt install -y python2 python-tk curl

RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py 
RUN python2 get-pip.py 
RUN pip install matplotlib==2.0.2 numpy 
RUN ln -s /usr/bin/python2 /usr/bin/python

RUN Rscript -e 'remotes::install_version("R.utils",upgrade="never", version = "2.10.1")'
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
RUN Rscript -e 'remotes::install_github("dcgerard/updog", ref = "eaa231c869e8d2e4c1fd11dc28e53e6270c9d5c0")'

RUN Rscript -e "install.packages('filelock')"
RUN Rscript -e "install.packages('vroom')"

RUN Rscript -e 'remotes::install_github("Cristianetaniguti/onemap")'

# Still in private repo
COPY ./onemapUTILS_0.3.1.tar.gz /opt/

RUN R -e "install.packages('/opt/onemapUTILS_0.3.1.tar.gz', type = 'source', repos=NULL)"


