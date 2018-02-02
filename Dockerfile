# Base image
FROM platform-docker.artifactory.niaid.nih.gov/shiny1-alpine3:latest
 
USER root
 
ENV JDK_VERSION="8"
ENV JAVA_HOME=/usr/lib/jvm/java-1.${JDK_VERSION}-openjdk

RUN apk update && apk upgrade && apk add --no-cache zip openjdk${JDK_VERSION}
 
RUN apk add --no-cache --virtual .build-dependencies make gcc R-dev g++ libxml2-dev && \
    mkdir -p /usr/share/doc/R/html && \
    R CMD javareconf && \
    R -e "install.packages('dplyr', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('leaflet', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('shinyjs', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('shinyBS', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('readr', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('stringi', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('reshape2', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('DT', repos='https://cran.rstudio.com/')" && \
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
    R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("org.Hs.eg.db", ask=FALSE); biocLite("org.Mm.eg.db", ask=FALSE);'  && \
    apk del .build-dependencies

# Make all files inside the directory 'app' available to the container
COPY app/ /srv/shiny-server/
RUN chown -R default /srv/shiny-server/
 
# Set or override the environmental variables
ENV LANG=en_US.UTF-8
ENV SHINY_APP_DATA_DIR  /srv/shiny-server/app/data/
ENV SHINY_APP_SCRIPT_DIR  /srv/shiny-server/app/Rscripts/
ENV SHINY_APP_INPUT_DIR  /srv/shiny-server/app/inputOutputs/TRIAGEinputFiles/
ENV SHINY_APP_OUTPUT_DIR  /srv/shiny-server/app/inputOutputs/TRIAGEoutputFiles/
 
# Override the previous user
USER default
