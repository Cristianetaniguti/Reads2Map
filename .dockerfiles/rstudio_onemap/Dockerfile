FROM cristaniguti/onemap_git:latest

RUN R -e "install.packages('R.utils', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('stringi', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('vcfR', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('doParallel', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('matrixStats', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('BiocManager', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('shiny', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('shinydashboard', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('tidyverse', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('largeList', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('ggpubr', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('config', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('irr', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('RColorBrewer', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('polyRAD', repos='http://brieger.esalq.usp.br/CRAN/');\
          BiocManager::install('pcaMethods');\
          BiocManager::install('VariantAnnotation');\
          devtools::install_github('tpbilton/GUSbase');\
          devtools::install_github('tpbilton/GUSMap')"

RUN R -e "devtools::install_github('Cristianetaniguti/genotyping4onemap', auth_token ='1bcc1d0962f96c9b4d9afc936b807ab2d0804d01')"

