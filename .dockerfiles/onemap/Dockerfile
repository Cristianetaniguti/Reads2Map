FROM r-base:3.6.0

RUN apt update \
   && apt-get install -y libcurl4-gnutls-dev libssl-dev libftgl2 libglu1-mesa-dev libftgl2 libfreetype6-dev  libgfortran5

RUN R -e "install.packages('updog', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('devtools', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('stringi', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('vcfR', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('reshape2', repos='http://brieger.esalq.usp.br/CRAN/' );\
          devtools::install_github('Cristianetaniguti/onemap_ht', auth_token ='1bcc1d0962f96c9b4d9afc936b807ab2d0804d01');\
          devtools::install_github('Cristianetaniguti/supermassa4onemap', auth_token ='1bcc1d0962f96c9b4d9afc936b807ab2d0804d01')"
