FROM ubuntu:18.04

RUN apt update \
    && apt install -y  build-essential wget \
    && wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.2/sratoolkit.2.10.2-ubuntu64.tar.gz \
    && tar xzf sratoolkit.2.10.2-ubuntu64.tar.gz \
    && cd sratoolkit.2.10.2-ubuntu64/bin \
    && ./vdb-config -i

WORKDIR /sratoolkit.2.10.2-ubuntu64/bin