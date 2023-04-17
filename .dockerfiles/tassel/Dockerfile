FROM ibmjava:jre

RUN apt update \
    && apt install -y  parallel vcftools

RUN conda install -y -c bioconda freebayes \
    && conda install -y -c bioconda vcflib

RUN wget https://bitbucket.org/tasseladmin/tassel-5-standalone/get/5f68583d0f56.zip \
    && unzip 5f68583d0f56.zip \
    && mv tasseladmin-tassel-5-standalone-5f68583d0f56/ /tassel
