##########################################
#      Ranking Hits From HuTNF
#       Based on Score and 
#     Network Interactions.
##########################################
# Samuel Katz
# 07/26/16
##########################################


library(data.table)                       #Load Libraries
library(dplyr)


#Set Working Directory

ScreenDir <- "~/Documents/Analysis/April_MoTNFtests/ImprovedAnnotationPlusNonKEGGaddition/"
SecScreenDir <- "~/Desktop/Secondary_Screens/LowestFourFracNeg/"
# TRIAGE output file goes here
AnalysisDir <- "~/Desktop/Sam_TRIAGE/Test_TRIAGErCode/TRIAGEoutputFiles/"
KEGGdir <- "~/Desktop/CARDcode/Rscripts/Resources/Pathways"
outDir <- "~/Documents/Analysis/Spliceosome/ScreenHits/"
TLR.dir <- "~/Documents/Analysis/GeneListsDB/"

#Get Edge, Node files, KEGG Files, TLR files, Screen Files
setwd(AnalysisDir)
EdgeInfo <- read.csv("HuTNF_D3NetworkEdgePropertiesHit.csv", stringsAsFactors = F)
NodeInfo <- read.csv("HuTNF_D3NetworkNodePropertiesHit.csv", stringsAsFactors = F)

setwd(KEGGdir)
KEGGhuman <- read.csv("KEGGHuman.csv", stringsAsFactors = F)

setwd(TLR.dir)
Ian.tlr.canon <- read.csv("TLR128.csv", stringsAsFactors = F)

setwd(AnalysisDir)
# Read in the TRIAGE outfile
ScreenScores <- read.csv("TRIAGEoutput_HuTNF_top5percCUTOFF_KEGG_hSTRINGppi.hi.csv", stringsAsFactors = F)
IAM_final_iteration <- colnames(ScreenScores)[(ncol(ScreenScores))-3]    #Define column from which to take "hit" data
InflectionPoint <- -4.5
# setwd(SecScreenDir)
# SecondaryScreen <- read.csv("lowestThreeAverageHuTNF_Entrez.csv", stringsAsFactors = FALSE)
# cutoff <- -1.32 #choose cutoff for Secondary Screen

#Get Matrix of genes for each pathway in KEGG
TLRpathway.genes.matrix <- matrix(toupper(Ian.tlr.canon$Human.Symbol[1:128]))

PROTpathway.genes <- filter(KEGGhuman, PathwayName == "Proteasome")
PROTpathway.genes.matrix <- matrix(PROTpathway.genes$GeneSymbol)

SPLICEpathway.genes <- filter(KEGGhuman, PathwayName == "Spliceosome")
SPLICEpathway.genes.matrix <- matrix(SPLICEpathway.genes$GeneSymbol)

#Add GeneSymbols to Edge dataframe
Edge.source <- merge(EdgeInfo, NodeInfo[, c("GeneMappingID", "GeneSymbol")], by.x = "source", by.y = "GeneMappingID", all.x = TRUE)
names(Edge.source)[names(Edge.source)=="GeneSymbol"] <- "source.ID"

Edge.target <- merge(Edge.source, NodeInfo[, c("GeneMappingID", "GeneSymbol")], by.x = "target", by.y = "GeneMappingID", all.x = TRUE)
names(Edge.target)[names(Edge.target)=="GeneSymbol"] <- "target.ID"

#Set column for pathway groupings
NodeInfo$Group <- "Novel"
NodeInfo$Group[NodeInfo$GeneSymbol %in% PROTpathway.genes.matrix] <- "Proteasome"
NodeInfo$Group[NodeInfo$GeneSymbol %in% SPLICEpathway.genes.matrix] <- "Splicesome"
NodeInfo$Group[NodeInfo$GeneSymbol %in% TLRpathway.genes.matrix] <- "TLRpathway"

#Add in column in Screen File for Hit y/n
ScreenScores <- ScreenScores %>% mutate(Hit_IAM = ifelse((Zscore < InflectionPoint & KEGGdb == "Absent" & CSAdes == "HitbyCSA" & get(IAM_final_iteration, envir = as.environment(ScreenScores)) != 1) | 
                                 get(IAM_final_iteration, envir = as.environment(ScreenScores)) == 1, 1, 0))

#Merge Zscore, Pathway, and Hit_IAM to NodeInfo
Scores_and_nodes <- merge(NodeInfo[, c("GeneMappingID", "GeneSymbol", "Group")], 
                          ScreenScores[, c("GeneSymbol", "EntrezID", "Zscore", "Pathway", "Hit_IAM")], 
                          by.x = "GeneSymbol", by.y = "GeneSymbol", all.x = T)

# #Add validation column to Secondary Screen file
# SecondaryScreen$Validation <- ifelse(SecondaryScreen$Replicate1 <= cutoff, 1, 0)
# 
# #Add column for presence in Secondary
# Scores_nodes_and_secondary <- Scores_and_nodes %>%
#   mutate(Sec.Screen = Scores_and_nodes$EntrezID %in% SecondaryScreen$EntrezID )
# 
# #convert values to zero and one.
# Scores_nodes_and_secondary$Sec.Screen[Scores_nodes_and_secondary$Sec.Screen == TRUE] <- 1
# Scores_nodes_and_secondary$Sec.Screen[Scores_nodes_and_secondary$Sec.Screen == FALSE] <- 0
# 
# #Merge Valdiation column
# Scores_nodes_and_secondary <- merge(Scores_nodes_and_secondary, 
#                                     SecondaryScreen[, c("EntrezID", "Validation")], 
#                                     by.x = "EntrezID", by.y = "EntrezID", all.x = T)
# 
# #Convert NA values to 0
# Scores_nodes_and_secondary$Validation[is.na(Scores_nodes_and_secondary$Validation)] <- 0

#Aggregate the edges to be summarised to each genemap ID
Edge_source_summary <- aggregate(target.ID ~ source.ID, data = Edge.target, paste, collapse = ", ")
Edge_target_summary <- aggregate(source.ID ~ target.ID, data = Edge.target, paste, collapse = ", ")


#Align column names and stack data frames
colnames(Edge_target_summary) <- c("source.ID", "target.ID")
Edge_summary_stacked <- rbind(Edge_source_summary, Edge_target_summary)

#Aggregate stacked data to get unique values
Edge_summary <- aggregate(target.ID ~ source.ID, data = Edge_summary_stacked, paste, collapse = ", ")

#Update Names
colnames(Edge_summary) <- c("GeneSymbol", "Ntwrk.all")

#Merge with scores and Secondary
Scores_nodes_secondary_and_edges <- merge(Scores_and_nodes, Edge_summary, by.x = "GeneSymbol", by.y = "GeneSymbol", all.x = T)

#Create column with counts for interacting genes
for (i in 1:length(Scores_nodes_secondary_and_edges$Ntwrk.all)) {
  Scores_nodes_secondary_and_edges$Allnet.count[i] <- ifelse(is.na(Scores_nodes_secondary_and_edges$Ntwrk.all[i]) == T, 0, 
                                                             length(unlist(strsplit(Scores_nodes_secondary_and_edges$Ntwrk.all[i], ", "))))
}

########################################### Create Matrices of the genes within each group that are also hits in the screen
#TLR Hit Genes 
TLRhits.genes <- filter(Scores_nodes_secondary_and_edges, Group == "TLRpathway", Hit_IAM == 1)
TLRhits.genes.matrix <- matrix(TLRhits.genes$GeneSymbol)

#Splicesome Hit Genes 
SPLICEhits.genes <- filter(Scores_nodes_secondary_and_edges, Group == "Splicesome", Hit_IAM == 1)
SPLICEhits.genes.matrix <- matrix(SPLICEhits.genes$GeneSymbol)

#TLR Hit Genes 
PROThits.genes <- filter(Scores_nodes_secondary_and_edges, Group == "Proteasome", Hit_IAM == 1)
PROThits.genes.matrix <- matrix(PROThits.genes$GeneSymbol)

########################################### Create Lists and Rankings by group------TLRpathway
#TLR Genes
#Create Blank
Scores_nodes_secondary_and_edges$Ntwrk.TLR <- NA
group.list <- TLRpathway.genes.matrix

#Extract genes from network interaction associted with genes from group
for (i in 1:length(group.list)) {
  temp.name <- group.list[i]
  out <- grep(paste0("\\b", temp.name, "\\b"), Scores_nodes_secondary_and_edges[["Ntwrk.all"]], value = F, fixed = F)
  Scores_nodes_secondary_and_edges$Ntwrk.TLR[out] <- paste(temp.name, Scores_nodes_secondary_and_edges$Ntwrk.TLR[out], sep = ", ")
}

#remove "NA" from rows with genes in them
Scores_nodes_secondary_and_edges$Ntwrk.TLR <- gsub(", NA", "", Scores_nodes_secondary_and_edges$Ntwrk.TLR)

#Create column with counts for interacting genes
for (i in 1:length(Scores_nodes_secondary_and_edges$Ntwrk.TLR)) {
  Scores_nodes_secondary_and_edges$TLRnet.count[i] <- ifelse(is.na(Scores_nodes_secondary_and_edges$Ntwrk.TLR[i]) == T, 0, 
                                                             length(unlist(strsplit(Scores_nodes_secondary_and_edges$Ntwrk.TLR[i], ", "))))
}

########################################### Create Lists and Rankings by group----Proteasome
#PROT Genes
#Create Blank
Scores_nodes_secondary_and_edges$Ntwrk.PROT <- NA
group.list <- PROTpathway.genes.matrix

#Extract genes from network interaction associted with genes from group
for (i in 1:length(group.list)) {
  temp.name <- group.list[i]
  out <- grep(paste0("\\b", temp.name, "\\b"), Scores_nodes_secondary_and_edges[["Ntwrk.all"]], value = F, fixed = F)
  Scores_nodes_secondary_and_edges$Ntwrk.PROT[out] <- paste(temp.name, Scores_nodes_secondary_and_edges$Ntwrk.PROT[out], sep = ", ")
}

#remove "NA" from rows with genes in them
Scores_nodes_secondary_and_edges$Ntwrk.PROT <- gsub(", NA", "", Scores_nodes_secondary_and_edges$Ntwrk.PROT)

#Create column with counts for interacting genes
for (i in 1:length(Scores_nodes_secondary_and_edges$Ntwrk.PROT)) {
  Scores_nodes_secondary_and_edges$PROTnet.count[i] <- ifelse(is.na(Scores_nodes_secondary_and_edges$Ntwrk.PROT[i]) == T, 0, 
                                                              length(unlist(strsplit(Scores_nodes_secondary_and_edges$Ntwrk.PROT[i], ", "))))
}

########################################### Create Lists and Rankings by group----Spliceseome
#SPLICE Genes
#Create Blank
Scores_nodes_secondary_and_edges$Ntwrk.SPLICE <- NA
group.list <- SPLICEpathway.genes.matrix

#Extract genes from network interaction associted with genes from group
for (i in 1:length(group.list)) {
  temp.name <- group.list[i]
  out <- grep(paste0("\\b", temp.name, "\\b"), Scores_nodes_secondary_and_edges[["Ntwrk.all"]], value = F, fixed = F)
  Scores_nodes_secondary_and_edges$Ntwrk.SPLICE[out] <- paste(temp.name, Scores_nodes_secondary_and_edges$Ntwrk.SPLICE[out], sep = ", ")
}

#remove "NA" from rows with genes in them
Scores_nodes_secondary_and_edges$Ntwrk.SPLICE <- gsub(", NA", "", Scores_nodes_secondary_and_edges$Ntwrk.SPLICE)

#Create column with counts for interacting genes
for (i in 1:length(Scores_nodes_secondary_and_edges$Ntwrk.SPLICE)) {
  Scores_nodes_secondary_and_edges$SPLICEnet.count[i] <- ifelse(is.na(Scores_nodes_secondary_and_edges$Ntwrk.SPLICE[i]) == T, 0, 
                                                                length(unlist(strsplit(Scores_nodes_secondary_and_edges$Ntwrk.SPLICE[i], ", "))))
}

##################################### Create list and rankings by group hits ##############################


############################# TLR Hits Genes

#Create Blank
Scores_nodes_secondary_and_edges$Ntwrk.TLR.hits <- NA
group.list <- TLRhits.genes.matrix

#Extract genes from network interaction associted with genes from group
for (i in 1:length(group.list)) {
  temp.name <- group.list[i]
  out <- grep(paste0("\\b", temp.name, "\\b"), Scores_nodes_secondary_and_edges[["Ntwrk.all"]], value = F, fixed = F)
  Scores_nodes_secondary_and_edges$Ntwrk.TLR.hits[out] <- paste(temp.name, Scores_nodes_secondary_and_edges$Ntwrk.TLR.hits[out], sep = ", ")
}

#remove "NA" from rows with genes in them
Scores_nodes_secondary_and_edges$Ntwrk.TLR.hits <- gsub(", NA", "", Scores_nodes_secondary_and_edges$Ntwrk.TLR)

#Create column with counts for interacting genes
for (i in 1:length(Scores_nodes_secondary_and_edges$Ntwrk.TLR.hits)) {
  Scores_nodes_secondary_and_edges$TLRhits.net.count[i] <- ifelse(is.na(Scores_nodes_secondary_and_edges$Ntwrk.TLR.hits[i]) == T, 0, 
                                                                  length(unlist(strsplit(Scores_nodes_secondary_and_edges$Ntwrk.TLR[i], ", "))))
}

############################# Proteasome Hits Genes

#Create Blank
Scores_nodes_secondary_and_edges$Ntwrk.PROT.hits <- NA
group.list <- PROThits.genes.matrix

#Extract genes from network interaction associted with genes from group
for (i in 1:length(group.list)) {
  temp.name <- group.list[i]
  out <- grep(paste0("\\b", temp.name, "\\b"), Scores_nodes_secondary_and_edges[["Ntwrk.all"]], value = F, fixed = F)
  Scores_nodes_secondary_and_edges$Ntwrk.PROT.hits[out] <- paste(temp.name, Scores_nodes_secondary_and_edges$Ntwrk.PROT.hits[out], sep = ", ")
}

#remove "NA" from rows with genes in them
Scores_nodes_secondary_and_edges$Ntwrk.PROT.hits <- gsub(", NA", "", Scores_nodes_secondary_and_edges$Ntwrk.PROT.hits)

#Create column with counts for interacting genes
for (i in 1:length(Scores_nodes_secondary_and_edges$Ntwrk.PROT.hits)) {
  Scores_nodes_secondary_and_edges$PROThits.net.count[i] <- ifelse(is.na(Scores_nodes_secondary_and_edges$Ntwrk.PROT.hits[i]) == T, 0, 
                                                                   length(unlist(strsplit(Scores_nodes_secondary_and_edges$Ntwrk.PROT.hits[i], ", "))))
}

############################# Splicesome Hits Genes

#Create Blank
Scores_nodes_secondary_and_edges$Ntwrk.SPLICE.hits <- NA
group.list <- SPLICEhits.genes.matrix

#Extract genes from network interaction associted with genes from group
for (i in 1:length(group.list)) {
  temp.name <- group.list[i]
  out <- grep(paste0("\\b", temp.name, "\\b"), Scores_nodes_secondary_and_edges[["Ntwrk.all"]], value = F, fixed = F)
  Scores_nodes_secondary_and_edges$Ntwrk.SPLICE.hits[out] <- paste(temp.name, Scores_nodes_secondary_and_edges$Ntwrk.SPLICE.hits[out], sep = ", ")
}

#remove "NA" from rows with genes in them
Scores_nodes_secondary_and_edges$Ntwrk.SPLICE.hits <- gsub(", NA", "", Scores_nodes_secondary_and_edges$Ntwrk.SPLICE.hits)

#Create column with counts for interacting genes
for (i in 1:length(Scores_nodes_secondary_and_edges$Ntwrk.SPLICE.hits)) {
  Scores_nodes_secondary_and_edges$SPLICEhits.net.count[i] <- ifelse(is.na(Scores_nodes_secondary_and_edges$Ntwrk.SPLICE.hits[i]) == T, 0, 
                                                                     length(unlist(strsplit(Scores_nodes_secondary_and_edges$Ntwrk.SPLICE.hits[i], ", "))))
}


############################## Network Hits Total
Scores_nodes_secondary_and_edges <- Scores_nodes_secondary_and_edges %>%
  mutate(Total_Path_Hits.net.count = TLRhits.net.count + SPLICEhits.net.count + PROThits.net.count)

##Calculate total network with other IAM hits
###write file
setwd(AnalysisDir)
write.csv(Scores_nodes_secondary_and_edges, "Ranking_HumanTNFScreen.csv")