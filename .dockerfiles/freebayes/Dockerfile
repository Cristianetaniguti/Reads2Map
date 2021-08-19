FROM conda/miniconda2

RUN apt update \
    && apt install -y  parallel vcftools

RUN conda install -y -c bioconda freebayes \
    && conda install -y -c bioconda vcflib
