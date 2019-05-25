FROM r-base:3.6.0

RUN apt update

RUN apt-get -y install libcurl4-gnutls-dev libssl-dev libftgl2 

RUN R -e "options(repos = \
   list(CRAN = 'http://brieger.esalq.usp.br/CRAN/')); \
   install.packages('updog');\
   install.packages('devtools')"

RUN apt -y install  libglu1-mesa-dev 

RUN R -e "devtools::install_github('Cristianetaniguti/onemap_ht', auth_token = '1bcc1d0962f96c9b4d9afc936b807ab2d0804d01')"

   

