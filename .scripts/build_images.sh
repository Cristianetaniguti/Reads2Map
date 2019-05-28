cd .dockerfiles/java-in-the-cloud
docker build -t java-in-the-cloud:v1 .
cd -

cd .dockerfiles/pIRS
docker build -t pirs:v1 .
cd -

cd .dockerfiles/onemap
docker build -t onemap:v1 .
cd -
