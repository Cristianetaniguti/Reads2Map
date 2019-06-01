cd .dockerfiles/java-in-the-cloud
docker build -t java-in-the-cloud:v1 .
cd -

cd .dockerfiles/onemap
docker build -t onemap:v1 .
cd -

cd .dockerfiles/pirs-ddrad-cutadapt
docker build -t pirs-ddrad-cutadapt:v1 .
cd -

cd .dockerfiles/bwa-samtools
docker build -t bwa-samtools:v1 .
cd -

cd .dockerfiles/gatk-picard
docker build -t gatk-picard:v1 .
cd -

cd .dockerfiles/freebayes
docker build -t freebayes:v1 .
cd -

cd .dockerfiles/stacks
docker build -t stacks:v1 .
cd -