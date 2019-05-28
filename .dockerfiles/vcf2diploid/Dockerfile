FROM ubuntu:18.04

RUN apt update
RUN apt-get update && \
    apt-get install -y build-essential default-jdk wget unzip && \
    wget https://github.com/abyzovlab/vcf2diploid/archive/master.zip && \
    unzip master.zip && \
    cd vcf2diploid-master/ && \
    make

WORKDIR vcf2diploid-master/