FROM rocker/rstudio:3.6.2

RUN apt update \
   && apt-get install -y libcurl4-gnutls-dev libssl-dev libftgl2 libglu1-mesa-dev libftgl2 libfreetype6-dev  libgfortran5 python-pip python-tk libxml2-dev libjpeg-dev libbz2-dev liblzma-dev \
   && pip install matplotlib==2.0.2 numpy 

RUN R -e "install.packages('updog', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('R.utils', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('devtools', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('stringi', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('vcfR', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('reshape2', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('ggplot2', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('doParallel', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('matrixStats', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('BiocManager', repos='http://brieger.esalq.usp.br/CRAN/'); \
          BiocManager::install('pcaMethods');\
          BiocManager::install('VariantAnnotation' );\
          install.packages('polyRAD', repos='http://brieger.esalq.usp.br/CRAN/');\
          devtools::install_github('tpbilton/GUSbase');\
          devtools::install_github('tpbilton/GUSMap')"

RUN R -e "devtools::install_github('Cristianetaniguti/onemap_ht', auth_token ='1bcc1d0962f96c9b4d9afc936b807ab2d0804d01', ref='parallel_map'  )"

RUN R -e "devtools::install_github('Cristianetaniguti/supermassa4onemap', auth_token ='1bcc1d0962f96c9b4d9afc936b807ab2d0804d01')"
