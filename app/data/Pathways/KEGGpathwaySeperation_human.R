##############################
### Obtain new KEGG databases
### Novemebr 19, 2020
### Samuel Katz
#############################

library("data.table")
library('dplyr')


#Get membership lists from KEGG
KEGG2020_pathwayDB <- fread('http://rest.kegg.jp/link/pathway/hsa', header = F)                         #Human gene-pathway membership
KEGGid_GeneID_human <- fread('http://rest.kegg.jp/conv/hsa/ncbi-geneid', header = F)                    #Human gene Kegg ID to NCBI gene ID map
KEGGpathwayIDs_human <- fread('http://rest.kegg.jp/list/pathway/hsa', header = F)                       #Human Kegg pathway to pathway ID conversion

##give column names to all tables.

colnames(KEGG2020_pathwayDB) <- c("KEGGgeneID", "KEGGpathID")
colnames(KEGGid_GeneID_human) <- c("NCBIgeneID", "KEGGgeneID")
colnames(KEGGpathwayIDs_human) <- c("KEGGpathID", "KEGGPathwayName")

#Add column with shaved off "ncbi-geneid:"

WordLenght <- nchar("ncbi-geneid:")
KEGGid_GeneID_human <- KEGGid_GeneID_human %>%
  mutate(EntrezID = substr(NCBIgeneID, WordLenght + 1, WordLenght + nchar(NCBIgeneID)))

#Add column with shaved off " - Homo sapiens (human):" From pathway names

WordLenght <- nchar(" - Homo sapiens (human)")
KEGGpathwayIDs_human <- KEGGpathwayIDs_human %>%
  mutate(PathwayName = substr(KEGGPathwayName, 1, nchar(KEGGPathwayName) - WordLenght))

# Merge file with Pathwaynames membership file
KEGG2020_pathwayDB <- merge(KEGG2020_pathwayDB, KEGGpathwayIDs_human[, c("KEGGpathID", "PathwayName")], all.x = T)

#Merge file with ENtrezIDs
KEGG2020_pathwayDB_EI <- merge(KEGG2020_pathwayDB, KEGGid_GeneID_human[, c("KEGGgeneID", "EntrezID")], by.x = "KEGGgeneID", by.y = "KEGGgeneID") 

#Get Genesymbol
library('org.Hs.eg.db')
Entrez_SYMBOL <- as.data.frame(org.Hs.egSYMBOL)

KEGG2020_pathwayDB_EI_GS <- merge(KEGG2020_pathwayDB_EI, Entrez_SYMBOL, by.x = "EntrezID", by.y = "gene_id", all.x = T)

colnames(KEGG2020_pathwayDB_EI_GS)[colnames(KEGG2020_pathwayDB_EI_GS) == "symbol"] <- "GeneSymbol"

#Create PathwayID column
WordLenght <- nchar("path:hsa")
KEGG2020_pathwayDB_EI_GS <- KEGG2020_pathwayDB_EI_GS %>%
  mutate(PathwayID = as.numeric(substr(KEGGpathID, WordLenght + 1, WordLenght + nchar(KEGGpathID))))






#Pull out non-disease pathways

KEGG2020_pathwayDB_EI_GS_nonDisease <- KEGG2020_pathwayDB_EI_GS %>%
  filter(PathwayID < 05000)

# KEGG2020_pathwayDB_EI_GS_HIV <- KEGG2020_pathwayDB_EI_GS %>%   #Create one for HIV analysis that has biological processes and HIV pathway
#   filter(PathwayName == "Human immunodeficiency virus 1 infection")
# KEGG2020_pathwayDB_EI_GS_HIV <- rbind(KEGG2020_pathwayDB_EI_GS_HIV, KEGG2020_pathwayDB_EI_GS_nonDisease)


KEGG2020_pathwayDB_EI_GS_Disease <- KEGG2020_pathwayDB_EI_GS %>%
  filter(PathwayID >= 05000)

#################


#######Write csv files
KEGG2020_Human_All <- KEGG2020_pathwayDB_EI_GS[, c("EntrezID", "GeneSymbol", "PathwayID", "PathwayName")]
setwd("~/OneDrive\ -\ National\ Institutes\ of\ Health/TRIAGE/PathwayAnalysis/KEGG2017/PathwayCSVs/201119")
write.csv(KEGG2020_Human_All, "KEGG2020_Human_All.csv")

KEGG2020_Human_BiologicalProcesses <- KEGG2020_pathwayDB_EI_GS_nonDisease[, c("EntrezID", "GeneSymbol", "PathwayID", "PathwayName")]
setwd("~/OneDrive\ -\ National\ Institutes\ of\ Health/TRIAGE/PathwayAnalysis/KEGG2017/PathwayCSVs/201119")
write.csv(KEGG2020_Human_BiologicalProcesses, "KEGG2020_Human_BiologicalProcesses.csv")

KEGG2020_Human_Disease <- KEGG2020_pathwayDB_EI_GS_Disease[, c("EntrezID", "GeneSymbol", "PathwayID", "PathwayName")]
setwd("~/OneDrive\ -\ National\ Institutes\ of\ Health/TRIAGE/PathwayAnalysis/KEGG2017/PathwayCSVs/201119")
write.csv(KEGG2020_Human_Disease, "KEGG2020_Human_Disease.csv")

# KEGG2020_HIV <- KEGG2020_pathwayDB_EI_GS_HIV[, c("EntrezID", "GeneSymbol", "PathwayID", "PathwayName")]
# setwd("~/Documents/TRIAGE/PathwayAnalysis/KEGG2017/PathwayCSVs/230119")
# write.csv(KEGG2020_HIV, "KEGG2020_HIVplusBioPro.csv")
# 



