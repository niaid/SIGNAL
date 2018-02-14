#############################################################
#           Changes to use StringDB and Entrez
#           in Network Analysis
#############################################################

#############################################################
#           Input parameters from the webpage
#############################################################
# /usr/bin/Rscript NetworkAnalysis.R "" "mmu" "NB-network" -1.5 "yes" "yes" "yes" "walktrap" "spring" 1000 5
          
#############################################################
#                   Load the Graphs
#############################################################
#rm(G,Graph)

G <- Selected_STRINGnetwork.igraph
message("Networks Loaded")
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
  message(Direction)
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
message(paste("Hit genes in network", length(hit.Genes.Network)))
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

