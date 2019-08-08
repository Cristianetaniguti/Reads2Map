FROM ubuntu:18.04

RUN apt-get update \
    && apt-get install -y build-essential git-core cmake zlib1g-dev libncurses-dev patch \
    && git clone https://github.com/genome/bam-readcount.git \
    && cd bam-readcount/ \
    && cmake . \
    && make
 
 WORKDIR /bam-readcount/bin/