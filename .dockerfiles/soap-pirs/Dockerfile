FROM continuumio/anaconda:latest

RUN apt-get update \   
    && apt-get install -y git build-essential libz-dev libboost-all-dev \
    && git clone --branch v2.0.2 https://github.com/galaxy001/pirs.git \
    && cd pirs \
    && make 

RUN conda install -c bioconda soapaligner \
    && conda install -c bioconda soapcoverage 

RUN git clone https://github.com/BGI-shenzhen/BamDeal.git \
    && cd BamDeal \
    && chmod 755 ./bin/BamDeal_Linux

RUN apt-get install -y samtools