## This Rscript is used to regenerate/update the human and mouse
## protein coding gene lists which are used to add the 'background' genes

library(data.table)
options(warn=-1)

##### Human Genes #####
## Complete list of all protein-encoding genes in human genome
## https://www.genenames.org/download/statistics-and-files/
## Website says 19202 genes, but only 19194 genes in the downloaded file!!!
## ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_types/gene_with_protein_product.txt
humanGenes <- fread('ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_types/gene_with_protein_product.txt')
human_tmp <- humanGenes[,c(2,19)]
names(human_tmp) <- c('GeneSymbol', 'EntrezID')
humanGenes = human_tmp

# Generate the lookup file for human protein coding genes
write.table(humanGenes, file="HGNC_genes_with_protein_product_EntrezID_geneSymbole_lookup.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

##### Mouse Genes #####
## Complete list of all genes in mouse genome
## 1. Get the protein coding genes only
## MGI Marker Coordinates (tab-delimited) - 2. Marker Type => 'protein coding gene'
## http://www.informatics.jax.org/downloads/reports/MGI_MRK_Coord.rpt
mouseGenes1 <- fread('http://www.informatics.jax.org/downloads/reports/MGI_MRK_Coord.rpt')
mouseProGenes <- mouseGenes1[mouseGenes1$'2. Marker Type' == 'protein coding gene', 4]
names(mouseProGenes) <- 'GeneSymbol'

## 2. Get the EntrezID, GeneSymbol 
mouseGenes2 <- fread('http://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt')
mouseGenes3 <- mouseGenes2[,c(3,6)]
names(mouseGenes3) <- c('GeneSymbol', 'EntrezID')

## 3. Get only the protein coding genes
mouseGenes4 <- mouseGenes3[mouseGenes3$GeneSymbol %in% mouseProGenes$GeneSymbol,]
mouseGenes = mouseGenes4

# Generate the lookup file for mouse protein coding genes
write.table(mouseGenes, file="Mouse_genes_with_protein_product_EntrezID_geneSymbol_lookup.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)