ComputeEnrichment <- function(pathway, hits, non.hits, file.name, siRNA.Score, iteration)
{
 
  hits <- intersect(hits,pathway$EntrezID)
  non.hits <- intersect(non.hits,pathway$EntrezID)
  pathway <- pathway[pathway$EntrezID %in% union(hits,non.hits),]
  unique.pathways <- unique(as.character(pathway$PathwayName))
  p.val <- pathway.genes.number <- hit.genes <- hit.gene.names <- rep(NA,length(unique.pathways))
  for(i in 1:length(unique.pathways))
  {
    pathway.genes <- pathway$EntrezID[which(as.character(pathway$PathwayName) == unique.pathways[i])]
    contingency <- matrix(NA,nrow = 2, ncol = 2)
    contingency[1,1] <- length(intersect(pathway.genes,hits)) # pathway.genes.hits
    contingency[1,2] <- length(intersect(pathway.genes,non.hits)) # pathway.genes.non.hits
    contingency[2,1] <- length(setdiff(hits,pathway.genes)) # non.pathway.hits 
    contingency[2,2] <- length(setdiff(non.hits,pathway.genes)) # non.pathway.non.hits
    p.val[i] <- fisher.test(contingency,alternative = "greater")$p.value
    pathway.genes.number[i] <- contingency[1,1] + contingency[1,2]
    hit.genes[i] <- contingency[1,1]
    hit.gene.names[i] <- paste(unique(pathway$GeneSymbol[match(intersect(pathway.genes,hits),pathway$EntrezID)]),collapse = ", ")
  }
  length(unique.pathways)
  p.val.FDR <- p.adjust(p.val,method = "BH")
  p.val.FWER <- p.adjust(p.val,method = "bonferroni")
  message('unique')

  results <- data.frame(Pathway = unique.pathways, 
                        pVal = round(p.val, digits = 3), 
                        pValFDR = round(p.val.FDR, digits = 3), 
                        pValBonferroni = round(p.val.FWER, digits = 3), 
                        Genes = pathway.genes.number, 
                        HitGenes = hit.genes, 
                        HitGeneNames = hit.gene.names)
  
  results <- results[with(results, order(results$pValBonferroni,results$pValFDR,results$pVal)),]
  write.csv(results,file = paste0(file.name,".Enrichment_", iteration, ".csv"), row.names = FALSE)
  
  
  sigPathways <- results$Pathway[which(results$pVal < 0.055)]
  sigPathwaysGenes <- unique(pathway$EntrezID[pathway$PathwayName %in% sigPathways])
  if(length(sigPathwaysGenes) > 0){
    nonSigPathwaysGenes <- setdiff(pathway$EntrezID, sigPathwaysGenes)
  } else {
    nonSigPathwaysGenes <- pathway$EntrezID
  }
  
  tempPathwayGenes <- matrix("Missing",nrow(siRNA.Score))
  message("sigPathways")
  message(sigPathways, "\n")  
  message(file.name)
  tempPathwayGenes[siRNA.Score$EntrezID %in% sigPathwaysGenes] <- "Yes"	
  tempPathwayGenes[siRNA.Score$EntrezID %in% nonSigPathwaysGenes] <- "No"
  if(length(grep("KEGG", file.name))>0){
    siRNA.Score$KEGG <- tempPathwayGenes
  } else if(length(grep("REACTOME", file.name))>0){
    siRNA.Score$REACTOME <- tempPathwayGenes
  } else if(length(grep("GO", file.name))>0){
    siRNA.Score$GO[intersect(which(siRNA.Score$GO == "Missing"),which(tempPathwayGenes == "No"))] <- "No"
    siRNA.Score$GO[siRNA.Score$EntrezID %in% sigPathwaysGenes] <- "Yes"	

  }
  siRNA.Score
}