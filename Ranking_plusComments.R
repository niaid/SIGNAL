###################################################
## Setting up Network Analysis 
## And ranking draft
###################################################
# Samuel Katz
# July 27, 2017
###################################################

require(edgebundleR)                        #Load Libraries
library('igraph')
library('data.table')
library('dplyr')

#############################################################
#Set directories

KEGGdir <- "~/Desktop/CARDcode/Rscripts/Resources/Pathways"      # This directory should countain a document with the memebrship lists of genes in pathways
PlotDir <- "~/Documents/Analysis/Simulating_Collar_Plot/"        # Where to place the plot you are going to create
Database.dir <- "~/Documents/Analysis/GeneListsDB/"              # Where to place additional membership lists, in this case this is where I placed  my own list of "TLR" genes
HitsDir <- "~/Documents/Analysis/April_MoTNFtests/ImprovedAnnotationPlusNonKEGGaddition/Testing/STRINGdb/" #Where to place the output of your TRIAGE analysis
CARDdirectory <- "~/Desktop/CARDcode/"                           # Where you placed your download from CARD (as well as the igraphs from STRING, in this case they are in a folder within it "./Rscripts/Resources/Network/")

#############################################################
                                                                 # I create matrices of the three pathways I want to look at. Currently this first step -choosing the pathways- is done manually. 
                                                                 # I chose TLR pathway because of the screen being a TLR screen, the spliceosome, and the proteasome because they were the highest scoring pathways
                                                                 # pathways that stood out to us.
##Getting Gene Hit Lists

#KEGG and pathwyay files                        
setwd(Database.dir)                                              #This was a list of "canonical" TLR genes that was provided to me (By Iain Fraser)
Ian.tlr.canon <- read.csv("TLR128.csv", stringsAsFactors = F)

setwd(KEGGdir)
KEGGhuman <- read.csv("KEGGHuman.csv", stringsAsFactors = F)

#Get Matrix of genes for each pathway in KEGG                    # Putting the list of gene EntrezID for each pathway of interest into a matrix
#Get TLR Canonical Genes

TLRpathway.genes.matrix <- matrix(na.omit(Ian.tlr.canon$Human.EntrezGene.ID))

PROTpathway.genes <- filter(KEGGhuman, PathwayName == "Proteasome")
PROTpathway.genes.matrix <- matrix(PROTpathway.genes$EntrezID)

SPLICEpathway.genes <- filter(KEGGhuman, PathwayName == "Spliceosome")
SPLICEpathway.genes.matrix <- matrix(SPLICEpathway.genes$EntrezID)

#Get IAM hits                                                    # Getting the TRIAGE output - name is hardcoded -I was working with Human TNF screen.
setwd(HitsDir)
HuTNFanno <- read.csv("IAMoutput_HuTNF_CSAincl_hSTRINGppi.hi.csv", stringsAsFactors = F)

#IAM hits and definging last iteration column name &  inflection point
IAM_final_iteration <- colnames(HuTNFanno)[(ncol(HuTNFanno))-3] # This column corresponds to the last network analysis step (expnasion) of the TRIAGE analysis.
InflectionPoint <- -4.5                                         # THis is the cutoff I used to "save" gene that aren't annotated in the pathway database (KEGG), the lowest score - since we were looking 
                                                                # at things with negative scores the "highest" was the one with the lowest zscore - was ~-6.5 so I used -4.5 as a cutoff.
#Get filtered IAMhits                                           # This is where I put together what is considered a "hit" by TRIAGE (IAM). Any gene that had a score of 1 in the last network analysis step
                                                                # and any gene that is in top two standard deviations that isn't annotated in KEGG
IAMhits <- filter(HuTNFanno, (Zscore < InflectionPoint & KEGGdb == "Absent" & CSAdes == "HitbyCSA" & get(IAM_final_iteration, envir = as.environment(HuTNFanno)) != 1) | 
                    get(IAM_final_iteration, envir = as.environment(HuTNFanno)) == 1)

IAMhits.matrix <- matrix(IAMhits$EntrezID)                      # Created a matrix of all the genes that are "hits" 


#Get filtered IAM file + Pathway genes                          # To do the complete networks, cerating a matrix of hit genes + plus genes that are in the pathways of interest
IAM <- filter(HuTNFanno, (Zscore < InflectionPoint & KEGGdb == "Absent" & CSAdes == "HitbyCSA" & get(IAM_final_iteration, envir = as.environment(HuTNFanno)) != 1) | 
                get(IAM_final_iteration, envir = as.environment(HuTNFanno)) == 1 | 
                EntrezID %in% TLRpathway.genes.matrix 
              #| EntrezID %in% PROTpathway.genes.matrix        # I commented out the protein pathway here and below but it can be commented back in.            
              | EntrezID %in% SPLICEpathway.genes.matrix)

#Get matrices for Pathways (Hits and non hits seperately)      # Now creating seperate matrices for each pathway containing only the genes of that pathway that are ALSO hits in the screen
#Hits
TLRhits <- filter(IAMhits, EntrezID %in% TLRpathway.genes.matrix)
TLRhits.matrix <- matrix(TLRhits$EntrezID)

# PROThits <- filter(IAMhits, EntrezID %in% PROTpathway.genes.matrix)
# PROThits.matrix <- matrix(PROThits$EntrezID)

SPLICEhits <- filter(IAMhits, EntrezID %in% SPLICEpathway.genes.matrix)
SPLICEhits.matrix <- matrix(SPLICEhits$EntrezID)

#Non hits                                                     # Now creating matrices of each pathways genes thate are not hits.
TLRnonhits.matrix <-  matrix(setdiff(TLRpathway.genes.matrix, IAMhits$EntrezID))

# PROTnonhits.matrix <- matrix(setdiff(PROTpathway.genes.matrix, IAMhits$EntrezID))

SPLICEnonhits.matrix <- matrix(setdiff(SPLICEpathway.genes.matrix, IAMhits$EntrezID))

#############################################################
#           Setting up Network Databases and input
#############################################################

organism <- "Human"

setwd(CARDdirectory)                                # Selecting which network to use, h/m human/mouse, hi/med hi confidence/medium confidence, ppi = Protein-Protein 
                                                    # interactions (as opposed to the others that combine other interaciton sources) 
if (tolower(organism) == "human") {
  networkTypes <- c("hSTRINGhi", "hSTRINGmed", "hSTRINGppi.hi", "hSTRINGppi.med")               
} else if (tolower(organism) == "mouse") {
  networkTypes <- c("mSTRINGhi", "mSTRINGmed", "mSTRINGppi.hi", "mSTRINGppi.med")
}
networkType <- networkTypes[c(3)] #**                          # Setting which network file to use.

use.only.commnected.components = "yes"
show.dendogram <- TRUE


siRNA.Score.Original <- IAM                                   #The gene input list, in this case all the genes (TRIAGE hits + pathways) that were put together earlier
#############################################################
#                   Load Library
#############################################################
library("igraph")
set.seed(123)
#############################################################  #Graphs loaded from ./Rscripts/Resources/Network/ it's written so that graphs can be added together cumilatively though this isn't used here
#                   Load the Graphs
#############################################################
if(tolower(organism) == "human") 
{  
  if("hSTRINGhi" %in% networkType)
  {
    load("./Rscripts/Resources/Network/igraph.string.hu.hiConf.Rdata")
    if(exists("G")) {
      G <- graph.union(igraph.string.hu.hiConf,G)
    } else {G <- igraph.string.hu.hiConf}
  }
  if("hSTRINGmed" %in% networkType)
  {
    load("./Rscripts/Resources/Network/igraph.string.hu.medConf.Rdata")
    if(exists("G")) {
      G <- graph.union(G,igraph.string.hu.medConf)
    } else {G <- igraph.string.hu.medConf}
  }
  if("hSTRINGppi.hi" %in% networkType)
  {
    load("./Rscripts/Resources/Network/igraph.stringPPI.hu.hiConf.Rdata")
    if(exists("G")) {
      G <- graph.union(igraph.stringPPI.hu.hiConf,G)
    } else {G <- igraph.stringPPI.hu.hiConf}
  }
  if("hSTRINGppi.med" %in% networkType)
  {
    load("./Rscripts/Resources/Network/igraph.stringPPI.hu.medConf.Rdata")
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
    load("./Rscripts/Resources/Network/igraph.string.mo.hiConf.Rdata")
    if(exists("G")) {G <- graph.union(igraph.string.mo.hiConf,G)}
    else {G <- igraph.string.mo.hiConf}
  }
  if("mSTRINGmed" %in% networkType)
  {
    cat('2')
    load("./Rscripts/Resources/Network/igraph.string.mo.medConf.Rdata")
    if(exists("G")) {G <- graph.union(igraph.string.mo.medConf,G)}
    else {G <- igraph.string.med.medConf}
  }
  if("mSTRINGppi.hi" %in% networkType)
  {
    cat('4')
    load("./Rscripts/Resources/Network/igraph.stringPPI.mo.hiConf.Rdata")
    if(exists("G")) {G <- graph.union(igraph.stringPPI.mo.hiConf,G)}
    else {G <- igraph.stringPPI.mo.hiConf}
  }
  if("mSTRINGppi.med" %in% networkType)
  {
    cat('2')
    load("./Rscripts/Resources/Network/igraph.stringPPI.mo.medConf.Rdata")
    if(exists("G")) {G <- graph.union(igraph.stringPPI.mo.medConf,G)}
    else {G <- igraph.stringPPI.mo.medConf}
  }
  
}
print("Networks Loaded")

#############################################################     # Using the methods from CARD (only swithced to EntrezID over GeneSymbol)
#            PageRank Algorithm to All Genes
#############################################################               
# AVERAGE DUPLICATED ROWS
siRNA.Score <- siRNA.Score.Original[which(!is.na(siRNA.Score.Original$EntrezID)),]    
OverallDegree <- degree(G)
Screen_Genes.for.network.analysis <- intersect(siRNA.Score$EntrezID,V(G)$name)
Graph <- induced.subgraph(G,Screen_Genes.for.network.analysis)

#############################################################
#            Select sub-graphs from hit genes
#############################################################
Subset_Genes.for.network.analysis <- Screen_Genes.for.network.analysis

SubGraph <- induced.subgraph(Graph,Subset_Genes.for.network.analysis)

if(tolower(use.only.commnected.components) == "yes"){
  SubGraph <- induced.subgraph(SubGraph,names(which(degree(SubGraph) > 0)))
}

#############################################################
#            Calculation of Node properties
#############################################################
Temp.Articulation <- V(SubGraph)$name[articulation.points(SubGraph)]
Articulation <- rep(0, length(V(SubGraph)$name))
Articulation[match(Temp.Articulation, (V(SubGraph)$name))] <- 1

siRNA.Score.Formatted <- siRNA.Score

gNames <- V(SubGraph)$name
###############################################################################
#                 Create Network for D3 Rendering for Hit Genes
###############################################################################
GraphEdgesHitNames <- get.data.frame(SubGraph, what = "edges")
GraphEdgesHitNames <- GraphEdgesHitNames[!duplicated(GraphEdgesHitNames), ]
source <- target <- rep(NA, nrow(GraphEdgesHitNames))
tempGenes <- union(GraphEdgesHitNames$from,GraphEdgesHitNames$to)
for(i in 1:length(tempGenes)){
  source[which(GraphEdgesHitNames$from == tempGenes[i])] <- i-1
  target[which(GraphEdgesHitNames$to == tempGenes[i])] <- i-1
}

GraphEdgesHitNumber <- data.frame(source,target)

GraphNodesHit <- data.frame(GeneMappingID = rep(0:(length(tempGenes)-1)), EntrezID = tempGenes)

###############################################################################
#                 Add back GeneSymbol and Hit Designation to output
###############################################################################
GraphNodesHit <- merge(GraphNodesHit, IAM[, c("EntrezID", "GeneSymbol", "ZSdes", IAM_final_iteration)], #Here my input file had a column "ZSdes" which indicated wether the gene made the zscore cutoff independently of TRIAGE, you can skip that part
                       by.x = "EntrezID", by.y = "EntrezID", all.x = T)

#############**************************************************################    # Now getting a data frame for the edges and a data frame for the nodes

EdgeInfo <- GraphEdgesHitNumber
NodeInfo <- GraphNodesHit


