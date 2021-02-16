##############################
### Obtain new KEGG databases
### June 29 2020
### Samuel Katz
#############################

library("data.table")
library('dplyr')


#Get membership lists from KEGG
KEGG2020_pathwayDB <- fread('http://rest.kegg.jp/link/pathway/mmu', header = F)                         #Mouse gene-pathway membership
KEGGid_GeneID_Mouse <- fread('http://rest.kegg.jp/conv/mmu/ncbi-geneid', header = F)                    #Mouse gene Kegg ID to NCBI gene ID map
KEGGpathwayIDs_Mouse <- fread('http://rest.kegg.jp/list/pathway/mmu', header = F)                       #Mouse Kegg pathway to pathway ID conversion

##give column names to all tables.

colnames(KEGG2020_pathwayDB) <- c("KEGGgeneID", "KEGGpathID")
colnames(KEGGid_GeneID_Mouse) <- c("NCBIgeneID", "KEGGgeneID")
colnames(KEGGpathwayIDs_Mouse) <- c("KEGGpathID", "KEGGPathwayName")

#Add column with shaved off "ncbi-geneid:"

WordLenght <- nchar("ncbi-geneid:")
KEGGid_GeneID_Mouse <- KEGGid_GeneID_Mouse %>%
  mutate(EntrezID = substr(NCBIgeneID, WordLenght + 1, WordLenght + nchar(NCBIgeneID)))

#Add column with shaved off " - Homo sapiens (Mouse):" From pathway names

WordLenght <- nchar(" - Homo sapiens (Mouse)")
KEGGpathwayIDs_Mouse <- KEGGpathwayIDs_Mouse %>%
  mutate(PathwayName = substr(KEGGPathwayName, 1, nchar(KEGGPathwayName) - WordLenght))

# Merge file with Pathwaynames membership file
KEGG2020_pathwayDB <- merge(KEGG2020_pathwayDB, KEGGpathwayIDs_Mouse[, c("KEGGpathID", "PathwayName")], all.x = T)

#Merge file with ENtrezIDs
KEGG2020_pathwayDB_EI <- merge(KEGG2020_pathwayDB, KEGGid_GeneID_Mouse[, c("KEGGgeneID", "EntrezID")], by.x = "KEGGgeneID", by.y = "KEGGgeneID") 

#Get Genesymbol
library('org.Mm.eg.db')
Entrez_SYMBOL <- as.data.frame(org.Mm.egSYMBOL)

KEGG2020_pathwayDB_EI_GS <- merge(KEGG2020_pathwayDB_EI, Entrez_SYMBOL, by.x = "EntrezID", by.y = "gene_id", all.x = T)

colnames(KEGG2020_pathwayDB_EI_GS)[colnames(KEGG2020_pathwayDB_EI_GS) == "symbol"] <- "GeneSymbol"

#Create PathwayID column
WordLenght <- nchar("path:mmu")
KEGG2020_pathwayDB_EI_GS <- KEGG2020_pathwayDB_EI_GS %>%
  mutate(PathwayID = as.numeric(substr(KEGGpathID, WordLenght + 1, WordLenght + nchar(KEGGpathID))))



#Pull out non-disease pathways

KEGG2020_pathwayDB_EI_GS_nonDisease <- KEGG2020_pathwayDB_EI_GS %>%
  filter(PathwayID < 05000)


KEGG2020_pathwayDB_EI_GS_Disease <- KEGG2020_pathwayDB_EI_GS %>%
  filter(PathwayID >= 05000)

#################
#######Write csv files
KEGG2020_Mouse_All <- KEGG2020_pathwayDB_EI_GS[, c("EntrezID", "GeneSymbol", "PathwayID", "PathwayName")]
setwd("~/OneDrive\ -\ National\ Institutes\ of\ Health/TRIAGE/PathwayAnalysis/KEGG2017/PathwayCSVs/201119")
write.csv(KEGG2020_Mouse_All, "KEGG2020_Mouse_All.csv")

KEGG2020_Mouse_BiologicalProcesses <- KEGG2020_pathwayDB_EI_GS_nonDisease[, c("EntrezID", "GeneSymbol", "PathwayID", "PathwayName")]
setwd("~/OneDrive\ -\ National\ Institutes\ of\ Health/TRIAGE/PathwayAnalysis/KEGG2017/PathwayCSVs/201119")
write.csv(KEGG2020_Mouse_BiologicalProcesses, "KEGG2020_Mouse_BiologicalProcesses.csv")

KEGG2020_Mouse_Disease <- KEGG2020_pathwayDB_EI_GS_Disease[, c("EntrezID", "GeneSymbol", "PathwayID", "PathwayName")]
setwd("~/OneDrive\ -\ National\ Institutes\ of\ Health/TRIAGE/PathwayAnalysis/KEGG2017/PathwayCSVs/201119")
write.csv(KEGG2020_Mouse_Disease, "KEGG2020_Mouse_Disease.csv")




