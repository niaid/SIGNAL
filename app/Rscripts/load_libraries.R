
library(shiny)
library(shinyjs)
library(shinyBS)
library(readr)
library(dplyr)
library(stringi)
library(DT)
library(data.table)
library(igraph)
library(edgebundleR)
library(shinyAce)
library(mailR)
library(rJava)
library(networkD3)
library(visNetwork)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(reshape2)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(crosstalk)
library(htmltools)
library(stringr)

L = c('shiny', 'shinyjs', 'shinyBS', 'readr', 'dplyr', 'stringi', 'DT', 'data.table',
      'igraph', 'edgebundleR', 'shinyAce', 'networkD3', 'visNetwork',
      'AnnotationDbi', 'reshape2', 'ggplot2', 'tidyr', 'gridExtra', 'crosstalk', 
      'htmltools', 'stringr', 'org.Hs.eg.db', 'org.Mm.eg.db', 'mailR', 'rJava',)


package.check <- lapply(L, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    if(x != 'org.Hs.eg.db' & x != 'org.Mm.eg.db'){
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
    else{
      if("BiocManager" %in% installed.packages()){
        BiocManager::install(x, version = "3.8")
      }
      else{
        install.packages("BiocManager")
        BiocManager::install(x, version = "3.8")
      }
    }
  }
})


load_libraries <- function(L){
  for(i in L){
    err <- F
    if(i != 'org.Hs.eg.db' & i != 'org.Mm.eg.db'){
      tryCatch({library(i, character.only = T)},
               error = function(err){err <<- T})
      ifelse(err, install.packages(i), library(i, character.only = T))
      }
    else{
      tryCatch({library(i, character.only = T)},
               error = function(err){err <<- T})
      if(err){
        
      }
      else{
        library(i, character.only = T)
      }
    }
  }
}

