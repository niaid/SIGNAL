# #Test variables
# screen.dataframe <- read.csv("/Users/sakatz/Documents/Analysis/HIV/TRIAGEHIV_May2019/TRIAGE/TRIAGEinput/Files/brass.triageinput.csv", stringsAsFactors = F)
# ID.column <- "EntrezID"
# criteria.column <- "normalized.triage"
# highconf.criteria <- 1
# midconf.criteria <- 0.5
# criteria.setting <- "equal"
# enrichment.dataframe <- read.csv("/Users/sakatz/Documents/Analysis/HIV/TRIAGEAnalyis/Rscripts/TRIAGEscripts/resources/KEGG2019_Human_BiologicalProcesses.csv", stringsAsFactors = F) 
# enrichment.title <- "PathwayName"
# stat.test <- "pVal"
# test.cutoff <- 0.055
# #- netowrk igraphs 
# load("~/Documents/Analysis/HIV/TRIAGEAnalyis/Rscripts/TRIAGEscripts/resources/String10.5.human.experimental.highConf.igraph.Rdata")
# load("~/Documents/Analysis/HIV/TRIAGEAnalyis/Rscripts/TRIAGEscripts/resources/String10.5.human.database.highConf.igraph.Rdata")
# load("~/Documents/Analysis/HIV/TRIAGEAnalyis/Rscripts/TRIAGEscripts/resources/String10.5.human.experimental.midConf.igraph.Rdata")
# load("~/Documents/Analysis/HIV/TRIAGEAnalyis/Rscripts/TRIAGEscripts/resources/String10.5.human.database.midConf.igraph.Rdata")
# network.igraph <- graph.union(String.human.experimental.highConf.igraph, String.human.database.highConf.igraph, String.human.experimental.midConf.igraph, String.human.database.midConf.igraph)  

####################################---
## TRIAGE AS FUNCTION
## (OR SET OF FUNCTIONS)
####################################---


########## ENRICHMENT 2 FUNCTION ----

###** Requirments **###
## screen.datafame: A dataframe of the screen 
## ID.column: A column within the screen.dataframe for the identifiers of the targets (EntrezID, GeneSymbol, etc.) 
## criteria.column: A column within the screen.dataframe of the criteria for being considered a hit 
## highconf.criteria: A criteria each target has to meet to be considered a "high confidence" hit.
## midconf.criteria: A criteria each target has to meet to be considered a "mid confidence" hit.
## criteria.setting: Whether you should be using "equal", "greater than or equal", or "less then or equal". Should be in the format of "equal", "greater", or "less"
## enrichment.dataframe: A dataframe to be used for pathway membership in the format of a column of IDs (should be same as ID column in screen.dataframe in ID type and column title) and a column of which group they are part of of. (each ID~group relationship should be in its own seperate row)
## enrichment.title: Name of the column with the names of the enrichment groups the targets are members of.
## stat.test: name of of the statistical test to be used for measuring enrichment confidence. Should be in format of either "pVal", "FDR", or "Bonferroni"
## test.cutoff: A numeric value which a less than the value in stat.test will be considered a significant enrichment.


###** Output **###
## a list with 2 dataframes
## 1: the screen.dataframe with an appended column listing each ID if it is part of a significantly enriched enrichment ("Yes") and was either high or medium confidence in the input, if it is not ("No"), or if it is missing from the enrichment.dataframe ("Missing")
## 2: a dataframe listing all the enrichmen groups with number of members, number of hits, ID of hits, p value, FDR, and Boneferroni of enrichment.

ENRICHMENT.function <- function(screen.dataframe, ID.column, criteria.column, highconf.criteria, midconf.criteria, criteria.setting, enrichment.dataframe, enrichment.title, stat.test, test.cutoff){

  ##--Assign Dataframes
  #Get dataframe of hits and assign temp column names
  hits.df <- screen.dataframe[, c(which(colnames(screen.dataframe) == ID.column), which(colnames(screen.dataframe) == criteria.column))]
  names(hits.df) <- c("ID.temp", "criteria.temp")

  #Get dataframe of enrichment and assign temp names
  enrich.df <- enrichment.dataframe[, c(which(colnames(enrichment.dataframe) == ID.column), which(colnames(enrichment.dataframe) == enrichment.title))]
  names(enrich.df) <- c("ID.temp", "enrich.temp")

  ##--Get high confidence, medium confidence hits and non hits matrix
  if (criteria.setting == "equal") {
    highconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp == highconf.criteria)])
    midconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp == midconf.criteria)])
    all.hits <- union(highconf.hits, midconf.hits) # a seperate matrix of high and low confidence hits to be used for graph subsetting later
    non.hits <- setdiff(as.matrix(hits.df$ID.temp), all.hits)
  } else if (criteria.setting == "greater") {
    highconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp >= highconf.criteria)])
    midconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp >= midconf.criteria && hits.df$criteria.temp < highconf.criteria)])
    all.hits <- union(highconf.hits, midconf.hits) # a seperate matrix of high and low confidence hits to be used for graph subsetting later
    non.hits <- setdiff(as.matrix(hits.df$ID.temp), all.hits)
  } else if (criteria.setting == "less") {
    highconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp <= highconf.criteria)])
    midconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp <= midconf.criteria && hits.df$criteria.temp > highconf.criteria)])
    all.hits <- union(highconf.hits, midconf.hits) # a seperate matrix of high and low confidence hits to be used for graph subsetting later
    non.hits <- setdiff(as.matrix(hits.df$ID.temp), all.hits)
  } else {
    highconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp == highconf.criteria)])
    midconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp == midconf.criteria)])
    all.hits <- union(highconf.hits, midconf.hits) # a seperate matrix of high and low confidence hits to be used for graph subsetting later
    non.hits <- setdiff(as.matrix(hits.df$ID.temp), all.hits)
    message("criteria.setting not properly defined. Using 'equals to' default")
  }

  ##-- subset hits and non.hits list to only those that intersect with what is in the enrichment dataframe
  subset.hits <- intersect(highconf.hits, enrich.df$ID.temp)
  subset.all.hits <- intersect(all.hits, enrich.df$ID.temp)
  subset.non.hits <- intersect(union(non.hits, midconf.hits), enrich.df$ID.temp)

  ##- subset enrichment list for enrichments that have members in the screen.dataframe list and get a list of unique names
  enrich.filter <- enrich.df[which(enrich.df$ID.temp %in% hits.df$ID.temp), ]
  enrich.unique <- unique(as.character(enrich.filter$enrich.temp))


  ##- add create columns and rows for all the pathways
  p.val <- enrich.group.members.number <- enrich.hit.number <- enrich.hit.IDs <- rep(NA,length(enrich.unique))


  ##-- populate the columns with the values for each enrichment group

  #- Fischer's Test (greater then)
  for(i in 1:length(enrich.unique))
  {
    enirchment.members <- enrich.filter$ID.temp[which(as.character(enrich.filter$enrich.temp) == enrich.unique[i])]
    contingency <- matrix(NA,nrow = 2, ncol = 2)
    contingency[1,1] <- length(intersect(enirchment.members, subset.hits)) # pathway.genes.hits
    contingency[1,2] <- length(intersect(enirchment.members, subset.non.hits)) # pathway.genes.non.hits
    contingency[2,1] <- length(setdiff(subset.hits, enirchment.members)) # non.pathway.hits
    contingency[2,2] <- length(setdiff(subset.non.hits, enirchment.members)) # non.pathway.non.hits
    p.val[i] <- fisher.test(contingency, alternative = "greater")$p.value
    enrich.group.members.number[i] <- contingency[1,1] + contingency[1,2]
    enrich.hit.number[i] <- contingency[1,1]
    enrich.hit.IDs[i] <- paste(unique(enrich.filter$ID.temp[match(intersect(enirchment.members,subset.hits), enrich.filter$ID.temp)]),collapse = ", ")
  }

  p.val.FDR <- p.adjust(p.val,method = "BH") #Correction for multiple testing
  p.val.FWER <- p.adjust(p.val,method = "bonferroni") #Bonferroni Correction


  #- print number of enrichment groups being calculated
  message(paste0("unique enrichment groups being measured: ", length(enrich.unique)))

  ##-- Create enrichment results dataframe

  enrichment.results <- data.frame(Enrichment = enrich.unique,
                                   pVal = p.val,
                                   pValFDR = p.val.FDR,
                                   pValBonferroni = p.val.FWER,
                                   EnrichmentMembers = enrich.group.members.number,
                                   EnrichmentHitNumber = enrich.hit.number,
                                   EnrichmentHitID = enrich.hit.IDs)

  enrichment.results <- enrichment.results[with(enrichment.results, order(enrichment.results$pValBonferroni,enrichment.results$pValFDR,enrichment.results$pVal)),] #order dataframe from lowest to highest value of statistical test

  ###--- Assign individual IDs in the screen.dataframe if they are hits based on the enrichment analysis and the provided cutoff

  ##-- Get list of enrichment that are siginificantly enriched for
  if (stat.test == "pVal") {
    sig.enrichment <- as.character(enrichment.results$Enrichment[which(enrichment.results$pVal < test.cutoff)])
  } else if (stat.test == "FDR") {
    sig.enrichment <- as.character(enrichment.results$Enrichment[which(enrichment.results$pValFDR < test.cutoff)])
  } else if (stat.test == "Bonferroni") {
    sig.enrichment <- as.character(enrichment.results$Enrichment[which(enrichment.results$pValBonferroni < test.cutoff)])
  } else {
    sig.enrichment <- as.character(enrichment.results$Enrichment[which(enrichment.results$pVal < test.cutoff)])
    message("Statisitical test not properly defined. Unsing pVal as default")
  }

  #- list of hit IDs in significantly enriched enrichments
  sig.enrich.hitIDs <- intersect(unique(enrich.filter$ID.temp[which(enrich.filter$enrich.temp %in% sig.enrichment)]), subset.all.hits)

  #- list of IDs that are not members of the significantly
  if(length(sig.enrich.hitIDs > 0)){
    nonhits.sig.enrich.IDs <- setdiff(enrich.filter$ID.temp, sig.enrich.hitIDs)
  } else {
    nonhits.sig.enrich.IDs <- enrich.filter$ID.temp
  }

  ##-- append column to screen.dataframe of whether the row ID is annotated in the enrichment.dataframe and if it is or isn't part of a significantly enriched group
  temp.enrich.IDs <- matrix("Missing", nrow(hits.df))
  temp.enrich.IDs[hits.df$ID.temp %in% sig.enrich.hitIDs] <- "Yes"
  temp.enrich.IDs[hits.df$ID.temp %in% nonhits.sig.enrich.IDs] <- "No"

  screen.dataframe$Enrichment.hit <- temp.enrich.IDs



  ###--- Define return objects
  enrichment.output <- list(screen.dataframe, enrichment.results)

  ####---- Return
  return(enrichment.output)

}



########## NETWORK FUNCTION ----

###** Requirments **###
## screen.datafame: A dataframe of the screen 
## ID.column: A column within the screen.dataframe for the identifiers of the targets (EntrezID, GeneSymbol, etc.). 
## criteria.column: A column within the screen.dataframe of the criteria for being considered a high confidence hit, medium confidence hit, and low confidence/non hit.
## highconf.criteria: A criteria each target has to meet to be considered a "high confidence" hit.
## midconf.criteria: A criteria each target has to meet to be considered a "mid confidence" hit.
## criteria.setting: Whether you should be using "equal", "greater than or equal", or "less then or equal". Should be in the format of "equal", "greater", or "less".
## network.igraph: an igraph of the network to be used for network analysis (network igraph must use the same ID type as screen.dataframe)


###** Output **###
## The input screen.dataframe with two appended columns
## Network.analysis: A column with information on whether it was an "InputHighConfidenceHit" or "NetworkAnalysisAdded", if neither than "0" is assigned.
## Network.hit: A column on wether based on network enrichment the row ID is a hit, with designations "Yes" and "No"

NETWORK.function <- function(screen.dataframe, ID.column, criteria.column, highconf.criteria, midconf.criteria, criteria.setting, network.igraph){

  ##-- necessary libraries
  library("igraph")

  ##--Assign Dataframes
  #Get dataframe of hits and assign temp column names
  hits.df <- screen.dataframe[, c(which(colnames(screen.dataframe) == ID.column), which(colnames(screen.dataframe) == criteria.column))]
  names(hits.df) <- c("ID.temp", "criteria.temp")

  ##--Get high confidence, medium confidence hits and non hits matrix
  if (criteria.setting == "equal") {
    highconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp == highconf.criteria)])
    midconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp == midconf.criteria)])
    all.hits <- union(highconf.hits, midconf.hits) # a seperate matrix of high and low confidence hits to be used for graph subsetting later
    non.hits <- setdiff(as.matrix(hits.df$ID.temp), all.hits)
  } else if (criteria.setting == "greater") {
    highconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp >= highconf.criteria)])
    midconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp >= midconf.criteria && hits.df$criteria.temp < highconf.criteria)])
    all.hits <- union(highconf.hits, midconf.hits) # a seperate matrix of high and low confidence hits to be used for graph subsetting later
    non.hits <- setdiff(as.matrix(hits.df$ID.temp), all.hits)
  } else if (criteria.setting == "less") {
    highconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp <= highconf.criteria)])
    midconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp <= midconf.criteria && hits.df$criteria.temp > highconf.criteria)])
    all.hits <- union(highconf.hits, midconf.hits) # a seperate matrix of high and low confidence hits to be used for graph subsetting later
    non.hits <- setdiff(as.matrix(hits.df$ID.temp), all.hits)
  } else {
    highconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp == highconf.criteria)])
    midconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp == midconf.criteria)])
    all.hits <- union(highconf.hits, midconf.hits) # a seperate matrix of high and low confidence hits to be used for graph subsetting later
    non.hits <- setdiff(as.matrix(hits.df$ID.temp), all.hits)
    message("criteria.setting not properly defined. Using 'equals to' default")
  }

  ###--- Format igraph
  set.seed(123)

  ##-- Select subgraph that matches list of IDs in screen  ----> This will be the geneal all encompasing network
  G <- upgrade_graph(network.igraph)
  OverallDegree <- degree(G)
  screen.IDs.for.network.analysis <- intersect(hits.df$ID.temp, V(G)$name)
  Graph <- induced.subgraph(G, screen.IDs.for.network.analysis)

  ##-- Select subgraph that matches list of IDs in high confidence AND mid confidence hits ----> This will be the network of all "hits"
  Subset.allhitIDs.for.network.analysis <- intersect(all.hits, V(G)$name)
  allhits.SubGraph <- induced.subgraph(G, Subset.allhitIDs.for.network.analysis)
  allhits.SubGraph <- induced.subgraph(allhits.SubGraph, names(which(degree(allhits.SubGraph) > 0)))
  allhits.SubGraph.IDs <- V(allhits.SubGraph)$name

  ##-- Removes edges that don't connect a high confidence hit.

  graph.edges.allhits.names <- get.data.frame(allhits.SubGraph, what = "edges")
  indices.to.remove <- intersect(which((graph.edges.allhits.names$from %in% highconf.hits) == FALSE),
                                 which((graph.edges.allhits.names$to %in% highconf.hits) == FALSE))

  if(length(indices.to.remove) > 0){
    graph.edges.allhits.names <- graph.edges.allhits.names[-indices.to.remove, ]
  }

  ##-- Get net list of added High confidence hits
  new.highconf.SubGraph <- graph.data.frame(graph.edges.allhits.names, directed = FALSE, vertices = NULL)

  ##--Combine this list with the original high confidence hits
  new.highconf.IDs <- union(V(new.highconf.SubGraph)$name, highconf.hits)

  ###-- Create columns to append to screen.dataframe
  hits.df$Network.analysis <- 0
  hits.df$Network.analysis[match(new.highconf.IDs, hits.df$ID.temp)] <- "NetworkAnalysisAdded"
  hits.df$Network.analysis[match(highconf.hits, hits.df$ID.temp)] <- "InputHighConfidenceHit"
  hits.df$Network.hit <- "No"
  hits.df$Network.hit[match(new.highconf.IDs, hits.df$ID.temp)] <- "Yes"

  ##-- append to the screen.dataframe input
  screen.dataframe.out <- screen.dataframe
  screen.dataframe.out$Network.analysis <- hits.df$Network.analysis
  screen.dataframe.out$Network.hit <- hits.df$Network.hit

  ###--- Return the appended data frame.
  return(screen.dataframe.out)

}



########## TRIAGE FUNCTION ----


###** Libraries **###
library('dplyr')
library('data.table')
library('igraph')


###** Requirments **###
## screen.datafame: A dataframe of the screen 
## ID.column: A column within the screen.dataframe for the identifiers of the targets (EntrezID, GeneSymbol, etc.) 
## criteria.column: A column within the screen.dataframe of the criteria for being considered a hit 
## highconf.criteria: A criteria each target has to meet to be considered a "high confidence" hit.
## midconf.criteria: A criteria each target has to meet to be considered a "mid confidence" hit.
## criteria.setting: Whether you should be using "equal", "greater than or equal", or "less than or equal". Should be in the format of "equal", "greater", or "less"
## enrichment.dataframe: A dataframe to be used for pathway membership in the format of a column of IDs (should be same as ID column in screen.dataframe in ID type and column title) and a column of which group they are part of of. (each ID~group relationship should be in its own seperate row)
## enrichment.title: Name of the column with the names of the enrichment groups the targets are members of.
## stat.test: name of of the statistical test to be used for measuring enrichment confidence. Should be in format of either "pVal", "FDR", or "Bonferroni"
## test.cutoff: A numeric value which a less than the value in stat.test will be considered a significant enrichment.
## network.igraph: an igraph of the network to be used for network analysis (network igraph must use the same ID type as screen.dataframe)


###** Output **###
## Output is a list of 3 dataframes
## [[1]] input dataframe plus 'TRIAGE.hit' column, 
## [[2]] dataframe of high confidence and medium confidence designation at each iteration, 
## [[3]] dataframe of final TRIAGE enrichments


TRIAGE.function <- function(screen.dataframe, ID.column, criteria.column, highconf.criteria, midconf.criteria, criteria.setting, enrichment.dataframe, enrichment.title, stat.test, test.cutoff, network.igraph) {
  
  ###--- Save original inputs
  input.screen.dataframe <- screen.dataframe
  input.ID.column <- ID.column
  input.criteria.column <- criteria.column
  input.highconf.criteria <- highconf.criteria
  input.midconf.criteria <- midconf.criteria
  input.criteria.setting <- criteria.setting
  input.enrichment.dataframe <- enrichment.dataframe
  input.enrichment.title <- enrichment.title
  input.stat.test <- stat.test
  input.test.cutoff <- test.cutoff
  input.network.igraph <- network.igraph
  
  
  
  ###--- Define input high confidence hits, medium confidence hits, and background
  #Get dataframe of hits and assign temp column names
  hits.df <- screen.dataframe[, c(which(colnames(screen.dataframe) == ID.column), which(colnames(screen.dataframe) == criteria.column))]
  names(hits.df) <- c("ID.temp", "criteria.temp")
  
  ##-- Get high confidence, medium confidence hits and non hits matrix of the inputs
  if (criteria.setting == "equal") {
    input.highconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp == highconf.criteria)])
    input.midconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp == midconf.criteria)])
    input.all.hits <- union(input.highconf.hits, input.midconf.hits) # a seperate matrix of high and low confidence hits to be used for graph subsetting later
    input.non.hits <- setdiff(as.matrix(hits.df$ID.temp), input.all.hits)
  } else if (criteria.setting == "greater") {
    input.highconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp >= highconf.criteria)])
    input.midconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp >= midconf.criteria && hits.df$criteria.temp < highconf.criteria)])
    input.all.hits <- union(input.highconf.hits, input.midconf.hits) # a seperate matrix of high and low confidence hits to be used for graph subsetting later
    input.non.hits <- setdiff(as.matrix(hits.df$ID.temp), input.all.hits)
  } else if (criteria.setting == "less") {
    input.highconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp <= highconf.criteria)])
    input.midconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp <= midconf.criteria && hits.df$criteria.temp > highconf.criteria)])
    input.all.hits <- union(input.highconf.hits, input.midconf.hits) # a seperate matrix of high and low confidence hits to be used for graph subsetting later
    input.non.hits <- setdiff(as.matrix(hits.df$ID.temp), input.all.hits)
  } else {
    input.highconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp == highconf.criteria)])
    input.midconf.hits <- as.matrix(hits.df$ID.temp[which(hits.df$criteria.temp == midconf.criteria)])
    input.all.hits <- union(highconf.hits, input.midconf.hits) # a seperate matrix of high and low confidence hits to be used for graph subsetting later
    input.non.hits <- setdiff(as.matrix(hits.df$ID.temp), input.all.hits)
    input.message("criteria.setting not properly defined. Using 'equals to' default")
  }
  
  ###--- Begin setting the first iteration variable to 1 
  iteration <- 1
  
  #- Create dataframe for appending the dataframe through each iteration step
  append.hits.df <- hits.df
  names(append.hits.df)[which(colnames(append.hits.df) == "ID.temp")] <- input.ID.column  ## Necesary to rename the ID column to original name so that it aligns with the ID.column of the enrichment dataframe
  
  ####---- Set up a counter for while loop till iteration converge
  counter <- TRUE
  
  while (counter == TRUE) {
    
    ###--- Enrichment step ----> Contracting the hits set
    
    enrichment.step <- ENRICHMENT.function(screen.dataframe, ID.column, criteria.column, highconf.criteria, midconf.criteria, criteria.setting, enrichment.dataframe, enrichment.title, stat.test, test.cutoff)
    
    ##-- create enrichment output dataframe
    #- Get datframe from enrichemnt step
    enrichment.output.df <- enrichment.step[[1]]
    
    ##-- Append Data frame to be used for Network function input
    append.hits.df <- data.frame(append.hits.df, temp = "", stringsAsFactors = FALSE)
    Enrichment.iteration.name <- paste0("ENRICH.iteration_", iteration)
    names(append.hits.df)[names(append.hits.df) == "temp"] <- Enrichment.iteration.name
    append.hits.df[[Enrichment.iteration.name]][enrichment.output.df$Enrichment.hit == "Yes" & (append.hits.df[[ID.column]] %in% input.all.hits)] <- "HighConf"
    append.hits.df[[Enrichment.iteration.name]][enrichment.output.df$Enrichment.hit != "Yes" & (append.hits.df[[ID.column]] %in% input.all.hits)] <- "MedConf"
    
    
    ###--- Network step ----> Expanding the hits set
    
    ##-- Set Network Function input parameters
    screen.dataframe <- append.hits.df
    ID.column <- ID.column
    criteria.column <- Enrichment.iteration.name
    highconf.criteria <- "HighConf"
    midconf.criteria <- "MedConf"
    criteria.setting <- "equal"
    network.igraph <- network.igraph
    
    
    ##-- Run NETWORK function
    
    network.output <- NETWORK.function(screen.dataframe, ID.column, criteria.column, highconf.criteria, midconf.criteria, criteria.setting, network.igraph)
    network.output.df <- network.output[[1]]
    
    ##-- Append Data frame to be used as Enrichment input or final output
    append.hits.df <- data.frame(append.hits.df, temp = "", stringsAsFactors = FALSE)
    Network.iteration.name <- paste0("NETWORK.iteration_", iteration)
    names(append.hits.df)[names(append.hits.df) == "temp"] <- Network.iteration.name
    append.hits.df[[Network.iteration.name]][network.output.df$Network.hit == "Yes" & (append.hits.df[[ID.column]] %in% input.all.hits)] <- "HighConf"
    append.hits.df[[Network.iteration.name]][network.output.df$Network.hit != "Yes" & (append.hits.df[[ID.column]] %in% input.all.hits)] <- "MedConf"
    
    ##-- Set Enrichment Function input parameters for next iteration
    screen.dataframe <- append.hits.df
    ID.column <- input.ID.column
    criteria.column <- Network.iteration.name
    highconf.criteria <- "HighConf"
    midconf.criteria <- "MedConf"
    criteria.setting <- "equal"
    enrichment.dataframe <- enrichment.dataframe 
    enrichment.title <-  enrichment.title
    stat.test <- stat.test
    test.cutoff <- test.cutoff
    
    ##-- Print message on completion of iteration
    message(paste("iteration ", iteration, " Complete"))
    
    ##-- Measure if there is an iterating pattern
    converge.sequence <- 0
    
    if (iteration >=3) {
      highconf.length <- c(length(append.hits.df[[ID.column]][append.hits.df[[Network.iteration.name]]== "HighConf"]))
      
      for ( t in 1:iteration){
        if (identical(append.hits.df[[Network.iteration.name]], append.hits.df[[paste0("NETWORK.iteration_", iteration-t)]])){
          converge.sequence <- t
          break
        } else {
          highconf.length <- c(highconf.length, length(append.hits.df[[ID.column]][append.hits.df[[paste0("NETWORK.iteration_", iteration-t)]] == "HighConf"]))
        }
      }
    }
    
    ##-- See if Iteration is converging on output 
    
    if((iteration != 1 && identical(append.hits.df[[Network.iteration.name]], append.hits.df[[paste0("NETWORK.iteration_", iteration-1)]])) 
       || (converge.sequence > 0 
           && (length(append.hits.df[[ID.column]][append.hits.df[[Network.iteration.name]]== "HighConf"]) == max(highconf.length)) 
       ))  {
      ##-- Set counter to false
      counter <- FALSE
    } else {
      ##-- Update itertion number
      iteration <- iteration + 1
    }
    
  }
  ####---- end of while loop, Print message on completion of TRIAGE iteration
  message(paste("TRIAGE iterations complete, number of iterations: ", iteration))
  
  
  
  
  ####---- Create TRIAGE output dataframes
  input.plus.triage.df <- input.screen.dataframe
  input.plus.triage.df$TRIAGE.hit <- network.output.df$Network.hit
  
  hits.by.iteration.df <- append.hits.df
  names(hits.by.iteration.df)[which(colnames(hits.by.iteration.df) == "criteria.temp")] <- input.criteria.column
  
  triage.enrichment.df <- enrichment.step[[2]]
  
  ##-- Create List of outputs
  triage.output <- list(input.plus.triage.df, hits.by.iteration.df, triage.enrichment.df)
  
  #- Print message on list content
  message("list contents: \n [[1]] input dataframe with 'TRIAGE.hit' column, \n [[2]] dataframe of high confidence and medium confidence designation at each iteration, \n [[3]] dataframe of final TRIAGE enrichments")
  
  
  ####---- return list
  return(triage.output)
}