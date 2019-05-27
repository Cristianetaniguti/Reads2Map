FROM ubuntu:18.04

RUN apt update
RUN apt -y install build-essential
RUN apt -y install default-jdk
RUN apt -y install wget
RUN apt -y install unzip

RUN wget https://github.com/abyzovlab/vcf2diploid/archive/master.zip
RUN unzip master.zip 
RUN cd vcf2diploid-master/ && make

WORKDIR vcf2diploid-master/