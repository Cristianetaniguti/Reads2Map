FROM ibmjava:jre

COPY ./PedigreeSim.jar /usr/jars/

RUN apt-get update \
    && apt-get install -y build-essential default-jdk wget unzip \
    && wget https://github.com/abyzovlab/vcf2diploid/archive/master.zip \
    && unzip master.zip \
    && cd vcf2diploid-master/ \
    && make \
    && mv vcf2diploid.jar /usr/jars/ 
 

RUN wget https://bitbucket.org/tasseladmin/tassel-5-standalone/get/5f68583d0f56.zip \
    && unzip 5f68583d0f56.zip \
    && mkdir /usr/tassel \
    && mv tasseladmin-tassel-5-standalone-5f68583d0f56/* /usr/tassel

RUN apt-get install -y samtools