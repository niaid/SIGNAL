# Base image
FROM platform-docker.artifactory.niaid.nih.gov/shiny1-alpine3:latest

# Set the user that will execute the RUN instruction(s)
USER root

# Install R package dependencies for the application. 
# The system packages required for building the R packages are installed but then removed at the end in order to reduce the container size
RUN apk update && \
    apk upgrade && \
    apk add zip && \
    apk add --no-cache --virtual .build-dependencies R-dev g++ libxml2-dev && \
    R -e "install.packages('dplyr', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('leaflet', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('DT', repos='https://cran.rstudio.com/')" && \    
    R -e "install.packages('shinyjs', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('shinyBS', repos='https://cran.rstudio.com/')" && \
    R -e "install.packages('readr', repos='https://cran.rstudio.com/')" && \ 
    R -e "install.packages('stringi', repos='https://cran.rstudio.com/')" && \ 
    R -e "install.packages('reshape2', repos='https://cran.rstudio.com/')" && \ 
    R -e "install.packages('data.table', repos='https://cran.rstudio.com/')" && \ 
    R -e "install.packages('edgebundleR', repos='https://cran.rstudio.com/')" && \ 
    R -e "install.packages('igraph', repos='https://cran.rstudio.com/')" && \ 
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
