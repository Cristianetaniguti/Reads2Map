FROM openjdk:7

COPY ./PedigreeSim.jar /usr/jars/

RUN apt-get update \
    && apt-get install -y build-essential default-jdk wget unzip \
    && wget https://github.com/abyzovlab/vcf2diploid/archive/master.zip \
    && unzip master.zip \
    && cd vcf2diploid-master/ \
    && make \
    && mv vcf2diploid.jar /usr/jars/ 
 
