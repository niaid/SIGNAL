FROM platform-docker.artifactory.niaid.nih.gov/shiny1-ubuntu18.04:latest 
 
USER root

RUN build_deps="r-base-dev openjdk-8-jdk libudunits2-dev libcairo2-dev libssl-dev libcurl4-openssl-dev" && \
    install_opts="-y --no-install-recommends" && \
    apt-get update && apt-get install $install_opts $build_deps && \
    R CMD javareconf && \
    R -e "install.packages('BiocManager')" && \
    R -e "BiocManager::install()" && \
    R -e "BiocManager::install('dplyr')" && \
    R -e "BiocManager::install('leaflet')" && \
    R -e "BiocManager::install('shinyjs')" && \
    R -e "BiocManager::install('shinyBS')" && \
    R -e "BiocManager::install('readr')" && \
    R -e "BiocManager::install('stringi')" && \
    R -e "BiocManager::install('stringr')" && \
    R -e "BiocManager::install('reshape2')" && \
    R -e "BiocManager::install('DT')" && \
    R -e "BiocManager::install('data.table')" && \
    R -e "BiocManager::install('edgebundleR')" && \
    R -e "BiocManager::install('igraph')" && \
    R -e "BiocManager::install('shinyAce')" && \
    R -e "BiocManager::install('rJava')" && \
    R -e "BiocManager::install('mailR')" && \
    R -e "BiocManager::install('networkD3')" && \
    R -e "BiocManager::install('visNetwork')" && \
    R -e "BiocManager::install('ggplot2')" && \
    R -e "BiocManager::install('tidyr')" && \
    R -e "BiocManager::install('gridExtra')" && \
    R -e "BiocManager::install('crosstalk')" && \
    R -e "BiocManager::install('htmltools')" && \
    R -e "BiocManager::install(c('org.Hs.eg.db', 'org.Mm.eg.db', 'AnnotationDbi'))" && \
    apt-get purge -y --auto-remove $build_deps && \
    apt-get install -y --no-install-recommends openjdk-8-jre && \ 
    rm -rf /var/lib/{apt,dpkg}
 
# Set or override the environmental variables
ENV SHINY_APP_DATA_DIR  /srv/shiny-server/app/data/
ENV SHINY_APP_SCRIPT_DIR  /srv/shiny-server/app/Rscripts/
ENV SHINY_APP_INPUT_DIR  /srv/shiny-server/app/inputOutputs/TRIAGEinputFiles/
ENV SHINY_APP_OUTPUT_DIR  /srv/shiny-server/app/inputOutputs/TRIAGEoutputFiles/
 
# Make all files inside the directory 'app' available to the container
COPY app/ /srv/shiny-server/
RUN chown -R default:default /srv/shiny-server

# Override the previous user
USER default
