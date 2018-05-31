FROM platform-docker.artifactory.niaid.nih.gov/shiny1-ubuntu18.04:latest 
 
USER root

RUN build_deps="r-base-dev openjdk-8-jdk libudunits2-dev libcairo2-dev libssl-dev libcurl4-openssl-dev" && \
    install_opts="-y --no-install-recommends" && \
    apt-get update && apt-get install $install_opts $build_deps && \
    R CMD javareconf && \
    R -e "install.packages('dplyr', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('leaflet', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('shinyjs', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('shinyBS', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('readr', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('stringi', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('reshape2', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('DT', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('data.table', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('edgebundleR', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('igraph', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('shinyAce', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('rJava', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('mailR', dep=T, repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('networkD3', dep=T, repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('visNetwork', dep=T, repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('ggplot2', dep=T, repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('tidyr', dep=T, repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('gridExtra', dep=T, repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('crosstalk', dep=T, repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('htmltools', dep=T, repos='https://cran.rstudio.com/')" && \
    R -e 'source("https://bioconductor.org/biocLite.R"); biocLite("org.Hs.eg.db", ask=FALSE); biocLite("org.Mm.eg.db", ask=FALSE); biocLite("AnnotationDbi", ask=FALSE);' && \
    apt-get purge -y --auto-remove $build_deps && \
    apt-get install -y --no-install-recommends openjdk-8-jre && \ 
    rm -rf /var/lib/{apt,dpkg}
 
# Set or override the environmental variables
ENV SHINY_APP_DATA_DIR  /srv/shiny-server/app/data/
ENV SHINY_APP_SCRIPT_DIR  /srv/shiny-server/app/Rscripts/
ENV SHINY_APP_INPUT_DIR  /srv/shiny-server/app/inputOutputs/TRIAGEinputFiles/
ENV SHINY_APP_OUTPUT_DIR  /srv/shiny-server/app/inputOutputs/TRIAGEoutputFiles/
 
# Override the previous user
USER default
# Make all files inside the directory 'app' available to the container
COPY app/ /srv/shiny-server/
