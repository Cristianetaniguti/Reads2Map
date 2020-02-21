FROM ubuntu:18.04

RUN apt update \
    && apt install -y git build-essential zlib1g-dev libbz2-dev liblzma-dev parallel python \
    && git clone --recursive --branch v1.3.0 https://github.com/ekg/freebayes.git \
    && cd freebayes \
    && make \
    && make install


RUN cd /freebayes/vcflib/ \
    && make

ENV PATH="/freebayes/vcflib/bin:${PATH}"
ENV PATH="/freebayes/scripts:${PATH}"
ENV PATH="/freebayes/vcflib/scripts:${PATH}"