##################################################################################################
## NAME: TRIAGE.R
## DESCRIPTION: Iterate microarray screen data through Pathway analysis and Network Analysis.
## REQUIREMENTS: 'Network_iteration_V2.R', 'Pathway_iteration.R'
## INPUTS: 
##      - Pathway Analysis: NormalizedData.csv 
## OUTPUTS:
##      - 'siRNA.Score'
##
## AUTHOR(S): Samuel Katz, Jian Song
## INSTITUTION: NIAID/NIH
## DATE LAST MODIFIED: 08/11/2017
#################################################################################################
scriptDir <- "~/TRIAGE/Rscripts/"
inputDir <-"~/TRIAGE/inputOutputs/TRIAGEinputFiles/"
outputDir <- "~/TRIAGE/inputOutputs/TRIAGEoutputFiles/"
dataDir <- "~/data/"
source(paste0(scriptDir, "pathway_iteration.R"))
Organism <- "Human"
network.type <- "hSTRINGppi.hi"

outDir <- "~/TRIAGE/inputOutputs/TRIAGEoutputFiles"
outputFileName <- paste0("TRIAGEoutput_HuTNF_CSAfdr_top5percCUTOFF_KEGG_", network.type, ".csv")

# 1) Seed Pathway Analysis
setwd(dataDir)


pathwayAnalysis <- function(pathway, input, cutoff){




pathway.types <- c("KEGG", "Reactome", "Gene_Ontology")
pathway.type <- pathway.types[1]
pathwayData <- read.csv(file = paste0(dataDir, "Pathways/", pathway.type, Organism, ".csv"))

seedName <- "TRIAGEinput_HuTNF_CSAfdr_5percCO.csv"

setwd(inputDir)
siRNA.Score <- read.csv(seedName, stringsAsFactors = F)
proxyScore <- "Replicate1" 
iteration <- 1
counter <- TRUE

while (counter == TRUE) {

  Hits <- siRNA.Score$EntrezID[siRNA.Score[[proxyScore]] == 1]
  nonHits <- setdiff(siRNA.Score$EntrezID, Hits)
  
  outPrefix <- paste0(pathway.type, iteration, sep = "_")  
  
  # 1) Contraction - [Pathway Analysis]
  siRNA.Score <- ComputeEnrichment(pathwayData, Hits, nonHits, outPrefix, siRNA.Score, iteration)
  siRNA.Score <- data.frame(siRNA.Score, temp = 0, stringsAsFactors = FALSE)
  kName1 <- paste0("KEGG.class.iteration", iteration)
  kName2 <- paste0("KEGG.", iteration)
  names(siRNA.Score)[names(siRNA.Score) == "temp"] <- kName1
  siRNA.Score[[kName1]][siRNA.Score$KEGG == "Yes" & siRNA.Score[[proxyScore]] > 0] <- 1
  siRNA.Score[[kName1]][siRNA.Score$KEGG != "Yes" & siRNA.Score[[proxyScore]] > 0] <- 0.5
  names(siRNA.Score)[names(siRNA.Score) == "KEGG"] <- kName2
  
  hit.Genes <- siRNA.Score$EntrezID[siRNA.Score[[kName1]] == 1]
  
  # 2) Expansion - [Network Analysis]
  system.time(
  source(paste0(scriptDir, "Network_iteration_V2.R"))
  )
  siRNA.Score <- data.frame(siRNA.Score, temp1 = "No", temp2 = siRNA.Score[[kName1]], stringsAsFactors = FALSE)
  nName1 <- paste0("Network.", iteration)
  nName2 <- paste0("Network.class.iteration", iteration)
  names(siRNA.Score)[names(siRNA.Score) == "temp1"] <- nName1
  names(siRNA.Score)[names(siRNA.Score) == "temp2"] <- nName2
  siRNA.Score[[nName1]][siRNA.Score$EntrezID %in% gNames2] <- "Yes"
  siRNA.Score[[nName2]][siRNA.Score$EntrezID %in% gNames2 & siRNA.Score[[kName1]] > 0] <- 1

  if(iteration != 1 && identical(siRNA.Score[[nName2]], siRNA.Score[[paste0("Network.class.iteration", iteration-1)]])) {
    dupCols <- (ncol(siRNA.Score)-1):ncol(siRNA.Score)
    siRNA.Score <- siRNA.Score[, -dupCols]
    counter <- FALSE
  }
  
  print(paste("iteration: ", iteration))
  proxyScore <- nName2
  iteration <- iteration + 1
}

### Append Enrichment Info --------------------------------------------
pval_threshold <- 0.05
enrichFileName <- paste0(outPrefix,".Enrichment_", iteration-1, ".csv")
pathEnrich <- read.csv(enrichFileName, stringsAsFactors = FALSE)

pathEnrich <- pathEnrich[pathEnrich$pVal < pval_threshold, ]

tempL <- strsplit(pathEnrich$HitGeneNames, split = ", ")
names(tempL) <- pathEnrich$Pathway

library(reshape2)
tempL <- lapply(seq(tempL), function(i) {
  m <- melt(tempL[i])
  names(m) <- c("GeneSymbol", paste0("Pathway", i))
  return(m)
})

tempDF <- Reduce(function(x, y) merge(x, y, by = "GeneSymbol", all = T), tempL)

samp <- tempDF[, 2:ncol(tempDF)]
pathVector <- sapply(seq(nrow(samp)), function(i) unlist(paste(samp[i, which(!is.na(samp[i, ]))], collapse = " ; ")))

pathDF <- data.frame(GeneSymbol = tempDF$GeneSymbol, Pathway = pathVector, stringsAsFactors = F)

out <- merge(siRNA.Score, pathDF, by = "GeneSymbol", all = T)

setwd(outDir)
write.csv(out, file = outputFileName, row.names = F)




#############################################################
#           Changes to use StringDB and Entrez
#           in Network Analysis
#############################################################



#############################################################
#           Input parameters from the webpage
#############################################################
# /usr/bin/Rscript NetworkAnalysis.R "" "mmu" "NB-network" -1.5 "yes" "yes" "yes" "walktrap" "spring" 1000 5
directory <- dataDir
organism <- Organism

setwd(directory)           
if (tolower(organism) == "human") {
  networkTypes <- c("hSTRINGhi", "hSTRINGmed", "hSTRINGppi.hi", "hSTRINGppi.med")               
} else if (tolower(organism) == "mouse") {
  networkTypes <- c("mSTRINGhi", "mSTRINGmed", "mSTRINGppi.hi", "mSTRINGppi.med")
}
networkType <- networkTypes[c(3)] #**

use.only.commnected.components = "yes"
show.dendogram <- TRUE
#fileName.siRNA.Score <-  "~/Desktop/Sam_CARD/NormalizedData_network.csv"

#siRNA.Score <- read.csv(file = fileName.siRNA.Score, header = TRUE); names(siRNA.Score)[2] <- "Score"

library("igraph")
set.seed(123)
#############################################################
#                   Load the Graphs
#############################################################
#rm(G,Graph)
if(tolower(organism) == "human") 
{  
  if("hSTRINGhi" %in% networkType)
  {
    load("~/data/Networks/igraph.string.hu.hiConf.Rdata")
    if(exists("G")) {
      G <- graph.union(igraph.string.hu.hiConf,G)
    } else {G <- igraph.string.hu.hiConf}
  }
  if("hSTRINGmed" %in% networkType)
  {
    load("~/data/Networks/igraph.string.hu.medConf.Rdata")
    if(exists("G")) {
      G <- graph.union(G,igraph.string.hu.medConf)
    } else {G <- igraph.string.hu.medConf}
  }
  if("hSTRINGppi.hi" %in% networkType)
  {
    load("~/data/Networks/igraph.stringPPI.hu.hiConf.Rdata")
    if(exists("G")) {
      G <- graph.union(igraph.stringPPI.hu.hiConf,G)
    } else {G <- igraph.stringPPI.hu.hiConf}
  }
  if("hSTRINGppi.med" %in% networkType)
  {
    load("~/data/Networks/igraph.stringPPI.hu.medConf.Rdata")
    if(exists("G")) {
      G <- graph.union(G,igraph.stringPPI.hu.medConf)
    } else {G <- igraph.stringPPI.hu.medConf}
  }

} else if(tolower(organism) == "mouse")
{
  cat('1')
  if("mSTRINGhi" %in% networkType)
  {
    cat('4')
    load("~/data/Networks/igraph.string.mo.hiConf.Rdata")
    if(exists("G")) {G <- graph.union(igraph.string.mo.hiConf,G)}
    else {G <- igraph.string.mo.hiConf}
  }
  if("mSTRINGmed" %in% networkType)
  {
    cat('2')
    load("~/data/Networks/igraph.string.mo.medConf.Rdata")
    if(exists("G")) {G <- graph.union(igraph.string.mo.medConf,G)}
    else {G <- igraph.string.med.medConf}
  }
  if("mSTRINGppi.hi" %in% networkType)
  {
    cat('4')
    load("~/data/Networks/igraph.stringPPI.mo.hiConf.Rdata")
    if(exists("G")) {G <- graph.union(igraph.stringPPI.mo.hiConf,G)}
    else {G <- igraph.stringPPI.mo.hiConf}
  }
  if("mSTRINGppi.med" %in% networkType)
  {
    cat('2')
    load("~/data/Networks/igraph.stringPPI.mo.medConf.Rdata")
    if(exists("G")) {G <- graph.union(igraph.stringPPI.mo.medConf,G)}
    else {G <- igraph.stringPPI.mo.medConf}
  }
  
}
print("Networks Loaded")
#############################################################
#            PageRank Algorithm to All Genes
#############################################################
# AVERAGE DUPLICATED ROWS
#siRNA.Score <- siRNA.Score[which(!is.na(siRNA.Score$Zscore)), ]     #commented out, as there is no column "zscore"

G <- upgrade_graph(G)
OverallDegree <- degree(G)
Screen_Genes.for.network.analysis <- intersect(siRNA.Score$EntrezID, V(G)$name)
Graph <- induced.subgraph(G, Screen_Genes.for.network.analysis)

########################################################################################################
#                               Identifying hit genes and apply off-target filters
########################################################################################################
findHitGenes <- function(x,Threshold, Direction, siRNAnumber = NULL){ # Function to find hit genes from a data.frame
  print(Direction)
  if(Direction == "Greater_than"){
    hit.Genes <- x$EntrezID[which(x[[proxyScore]] > Threshold)]   # Here gene symvols are used but in Network analysis gene IDs are used
  } else if(Direction == "Less_than"){
    hit.Genes <- x$EntrezID[which(x[[proxyScore]] < Threshold)]
  } else if(Direction == "Both"){
    hit.Genes <- x$EntrezID[union(which(x[[proxyScore]] > Threshold), which(x[[proxyScore]] < -1*Threshold))]
  }
  if(!is.null(siRNAnumber)){
    tempTable <- table(hit.Genes)
    hit.Genes <- names(tempTable[which(tempTable >= siRNAnumber)])
  }
  return(hit.Genes)
}
#############################################################
#            Select sub-graphs from hit genes
#############################################################
Subset_Genes.for.network.analysis <- intersect(hit.Genes, V(G)$name)
SubGraph <- induced.subgraph(G, Subset_Genes.for.network.analysis)
if(tolower(use.only.commnected.components) == "yes"){
  SubGraph <- induced.subgraph(SubGraph, names(which(degree(SubGraph) > 0)))
}
#############################################################
#            Calculation of Node properties
#############################################################
Temp.Articulation <- V(SubGraph)$name[articulation.points(SubGraph)]
Articulation <- rep(0, length(V(SubGraph)$name))
Articulation[match(Temp.Articulation, (V(SubGraph)$name))] <- 1

siRNA.Score.Formatted <- siRNA.Score

gNames <- V(SubGraph)$name


#############################################################
#                   Indirect Connection
#############################################################
NetworkConnectivityThreshold <- 1
NetworkScoreThreshold <- 0.4

tempPeripheralGenes <- PeripheralGenes <- NULL
hit.Genes.Network <- intersect(hit.Genes, V(Graph)$name)
print(paste("Hit genes in network", length(hit.Genes.Network)))
for(i in 1:length(hit.Genes.Network)){
  
  tempPeripheralGenes <- V(Graph)[unlist(neighborhood(Graph, 1, nodes=hit.Genes.Network[i], mode="all"))]$name
  if(length(tempPeripheralGenes) > 0) {
    tempPeripheralGenes <- setdiff(unique(tempPeripheralGenes), hit.Genes.Network)
  }
  PeripheralGenes <- c(tempPeripheralGenes, PeripheralGenes)
}
PeripheralGenesFrequency <- table(PeripheralGenes)
PeripheralGenesSelected <- names(PeripheralGenesFrequency[PeripheralGenesFrequency >= NetworkConnectivityThreshold])

tempGenes <- findHitGenes(siRNA.Score.Formatted, NetworkScoreThreshold, "Greater_than")
if(length(tempGenes) > 0) {
  PeripheralGenesSelected <- intersect(PeripheralGenesSelected, tempGenes)
}

PeripheralHitGenes <- union(hit.Genes.Network, PeripheralGenesSelected)
###############################################################################
#               Create subgraph and remove edges between non-hits 
###############################################################################
Subset_Genes.for.network.analysis <- intersect(PeripheralHitGenes,V(G)$name)
SubGraph <- induced.subgraph(G,Subset_Genes.for.network.analysis)
if(tolower(use.only.commnected.components) == "yes"){
  SubGraph <- induced.subgraph(SubGraph,names(which(degree(SubGraph) > 0)))
}
GraphEdgesHitPeripheralNames <- get.data.frame(SubGraph, what = "edges")
indRemove <- intersect(which((GraphEdgesHitPeripheralNames$from %in% gNames) == FALSE), 
                       which((GraphEdgesHitPeripheralNames$to %in% gNames) == FALSE))

if(length(indRemove) > 0){
  GraphEdgesHitPeripheralNames <- GraphEdgesHitPeripheralNames[-indRemove, ]
}
SubGraph <- graph.data.frame(GraphEdgesHitPeripheralNames, directed = FALSE, vertices = NULL)

#############################################################
#                  Create Subgraph with New Genes
#############################################################
Temp.Articulation <- V(SubGraph)$name[articulation.points(SubGraph)]
Articulation <- rep(0, length(V(SubGraph)$name))
Articulation[match(Temp.Articulation,(V(SubGraph)$name))] <- 1

gNames2 <- V(SubGraph)$name

ScreenHit <- rep("No", length(gNames2))
ScreenHit[match(gNames, gNames2)] <- "Yes"
indPeripheralGenes <- match(PeripheralGenesSelected, gNames2)
selected.gNames2 <- gNames2[indPeripheralGenes]


