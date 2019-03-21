### Changed for new Network, Feb 13, 2018, SK
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

#G <- Selected_STRINGnetwork.igraph
#############################################################
#            PageRank Algorithm to All Genes
#############################################################
# AVERAGE DUPLICATED ROWS
#siRNA.Score <- siRNA.Score[which(!is.na(siRNA.Score$Zscore)), ]     #commented out, as there is no column "zscore"
library("igraph")
set.seed(123)

OverallDegree <- degree(G)
Screen_Genes.for.network.analysis <- intersect(siRNA.Score$EntrezID, V(G)$name)
Graph <- induced.subgraph(G, Screen_Genes.for.network.analysis)

########################################################################################################
#                               Identifying hit genes and apply off-target filters
########################################################################################################
findHitGenes <- function(x, siRNAnumber = NULL){ # Function to find hit genes from a data.frame
  hit.Genes <- x$EntrezID[which(x[[proxyScore]] %in% c("MedConf", "HighConf"))]   # Here gene symvols are used but in Network analysis gene IDs are used
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
# NetworkConnectivityThreshold <- 1
# NetworkScoreThreshold <- 0.4

# NetworkConnectivityThreshold <- cutoffHigh
# NetworkScoreThreshold <- cutoffMed

tempPeripheralGenes <- PeripheralGenes <- NULL
hit.Genes.Network <- intersect(hit.Genes, V(Graph)$name)
message(paste("Hit genes in network", length(hit.Genes.Network)))

if(length(hit.Genes.Network) <= 0){
  # message("No hit gene found! Please change your cutoff values to ensure genes are selected.")
  # showNotification("No hit gene found! Please change your cutoff values to ensure genes are selected.", duration = NULL,
  #                  action = a(href = "javascript:location.reload();", "Reload page")
  # )
  showModal(modalDialog(
    title=HTML("<h3><font color=#ff0000>Error with cutoff values!</font></h3>"),
    HTML("No hit gene found! Please try changing your Interaction Sources or Confidence to ensure a set of genes will be selected!<br>
         Session will restart."),
    easyClose = TRUE
  ))
  Sys.sleep(5)
  session$reload()
} else{
  for(i in 1:length(hit.Genes.Network)){
    
    tempPeripheralGenes <- V(Graph)[unlist(neighborhood(Graph, 1, nodes=hit.Genes.Network[i], mode="all"))]$name
    if(length(tempPeripheralGenes) > 0) {
      tempPeripheralGenes <- setdiff(unique(tempPeripheralGenes), hit.Genes.Network)
    }
    PeripheralGenes <- c(tempPeripheralGenes, PeripheralGenes)
  }
  PeripheralGenesFrequency <- table(PeripheralGenes)
  PeripheralGenesSelected <- names(PeripheralGenesFrequency[PeripheralGenesFrequency >= 1])
  
  tempGenes <- findHitGenes(siRNA.Score.Formatted)
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

try(gNames2 <- V(SubGraph)$name, silent = FALSE)

ScreenHit <- rep("No", length(gNames2))
ScreenHit[match(gNames, gNames2)] <- "Yes"
indPeripheralGenes <- match(PeripheralGenesSelected, gNames2)
selected.gNames2 <- gNames2[indPeripheralGenes]
}

