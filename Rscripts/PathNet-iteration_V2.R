#############################################################
#           Changes to use StringDB and Entrez
#           in Network Analysis     #09/12/16
#############################################################


## NAME: PathNet-iteration.R
## DESCRIPTION: Iterate microarray screen data through Pathway analysis and Network Analysis.
## REQUIREMENTS: 'Network_iteration_V2.R', 'Pathway_iteration.R'
## INPUTS: 
##      - Pathway Analysis: NormalizedData.csv 
## OUTPUTS:
##      - 'siRNA.Score'
##
## AUTHOR(S): Bhaskar Dutta, Jason Guo, Samuel Katz, Nicolas Lounsbury
## INSTITUTION: NIAID/NIH
## DATE LAST MODIFIED: 02/09/2016
#################################################################################################


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
