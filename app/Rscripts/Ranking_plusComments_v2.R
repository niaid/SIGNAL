###################################################
## Setting up Network Analysis
## And ranking draft
###################################################
# Samuel Katz
# July 27, 2017
# modified by Jian Song
###################################################

# require('edgebundleR')                        #Load Libraries
# library('igraph')
# library('data.table')
# library('dplyr')

#selectedRows <- c(1,2,3)  
#Generate_NetworkGraph(selectedRows)

Generate_NetworkGraph <- function(selectedRows, organism){
  #############################################################
  #Set directories
  
  #TRIAGE.input <- "~/TRIAGE/app/data"
  TRIAGE.input <- dataDir

  #TRIAGE.output <- "~/TRIAGE/app/InputOutputs/TRIAGEoutputFiles"
  TRIAGE.output <- outputDir
  
  #Get pathway file
  pathwayFile <- pathwayData      
  
  
  
  #Get Matrix of genes for each pathway in KEGG                    # Putting the list of gene EntrezID for each pathway of interest into a matrix
  
  
  
  path1_name <<- sigPathways$Pathway[as.numeric(selectedRows[1])]     
  message(path1_name)
  path1_pathway.genes <- filter(pathwayFile, PathwayName == path1_name)
  path1_pathway.genes.matrix <- matrix(path1_pathway.genes$EntrezID)
  
  
  path2_name <<- sigPathways$Pathway[as.numeric(selectedRows[2])]   
  message(path2_name)
  path2_pathway.genes <- filter(pathwayFile, PathwayName == path2_name)
  path2_pathway.genes.matrix <- matrix(path2_pathway.genes$EntrezID)
  
  
  path3_name <<- sigPathways$Pathway[as.numeric(selectedRows[3])]    
  message(path3_name)
  path3_pathway.genes <- filter(pathwayFile, PathwayName == path3_name)
  path3_pathway.genes.matrix <- matrix(path3_pathway.genes$EntrezID)
  
  #Get TRIAGE Output                                             #Ultimately to be replaces with master data frame
  #TRIAGE.output.df <- read.csv(file = TRIAGE.output, stringsAsFactors = F)
  TRIAGE.output.df <- read.csv(paste0(outputDir, userDir, '/', outputFileName), stringsAsFactors = F)  #COMMMENT BACK IN FOR SHINY

  #IAM hits and definging last iteration column name &  inflection point
  TRIAGE_final_iteration <- colnames(TRIAGE.output.df)[(ncol(TRIAGE.output.df))-1] # This column corresponds to the last network analysis step (expnasion) of the TRIAGE analysis.
  
  #Get filtered TRIAGEhits                                           # This is where I put together what is considered a "hit" by TRIAGE (IAM). Any gene that had a score of 1 in the last network analysis step
  # and any gene that is in top two standard deviations that isn't annotated in KEGG
  TRIAGEhits <- filter(TRIAGE.output.df, get(TRIAGE_final_iteration, envir = as.environment(TRIAGE.output.df)) == 1)
  
  TRIAGEhits.matrix <- matrix(TRIAGEhits$EntrezID)                      # Created a matrix of all the genes that are "hits" 
  
  
  #Get matrices for Pathways (Hits and non hits seperately)      # Now creating seperate matrices for each pathway containing only the genes of that pathway that are ALSO hits in the screen
  #Hits
  path1_hits <- filter(TRIAGEhits, EntrezID %in% path1_pathway.genes.matrix)
  path1_hits.matrix <- matrix(path1_hits$EntrezID)
  
  path2_hits <- filter(TRIAGEhits, EntrezID %in% path2_pathway.genes.matrix)
  path2_hits.matrix <- matrix(path2_hits$EntrezID)
  
  path3_hits <- filter(TRIAGEhits, EntrezID %in% path3_pathway.genes.matrix)
  path3_hits.matrix <- matrix(path3_hits$EntrezID)
  
  #Non hits                                                     # Now creating matrices of each pathways genes thate are not hits.
  path1_nonhits.matrix <-  matrix(setdiff(path1_pathway.genes.matrix, TRIAGEhits$EntrezID))
  
  path2_nonhits.matrix <- matrix(setdiff(path2_pathway.genes.matrix, TRIAGEhits$EntrezID))
  
  path3_nonhits.matrix <- matrix(setdiff(path3_pathway.genes.matrix, TRIAGEhits$EntrezID))
  
  #############################################################
  #           Setting up Network Databases and input
  #############################################################
  
  
  use.only.commnected.components = "yes"
  show.dendogram <- TRUE
  
  
  #############################################################
  #                   Load Library
  #############################################################
  library("igraph")
  set.seed(123)
  #############################################################  #Graphs loaded from ~/TRIAGE/app/Rscripts/Resources/Network/ it's written so that graphs can be added together cumilatively though this isn't used here
  #                   Load the Graphs
  #############################################################
  G <- Selected_STRINGnetwork.igraph
  
  
  #############################################################     # Using the methods from CARD (only swithced to EntrezID over GeneSymbol)
  #            PageRank Algorithm to All Genes
  #############################################################               
  # AVERAGE DUPLICATED ROWS
  siRNA.Score <- TRIAGEhits[which(!is.na(TRIAGEhits$EntrezID)),]    
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
  GraphNodesHit <- merge(GraphNodesHit, TRIAGEhits[, c("EntrezID", "GeneSymbol", TRIAGE_final_iteration)], 
                         by.x = "EntrezID", by.y = "EntrezID", all.x = T)
  
  #############**************************************************################    # Now getting a data frame for the edges and a data frame for the nodes
  
  EdgeInfo <- GraphEdgesHitNumber
  NodeInfo <- GraphNodesHit
  
  
  #######################################                                            # redoing the pathway gene matrices to be of the gene symbol instead of the Entrez ID
  #Get Matrix of genes for each pathway in KEGG
  path1_pathway.genes.matrix <- matrix(path1_pathway.genes$GeneSymbol)
  
  path2_pathway.genes.matrix <- matrix(path2_pathway.genes$GeneSymbol)
  
  path3_pathway.genes.matrix <- matrix(path3_pathway.genes$GeneSymbol)
  
  
  #Add GeneSymbols to Edge dataframe                                              
  Edge.source <- merge(EdgeInfo, NodeInfo[, c("GeneMappingID", "GeneSymbol")], by.x = "source", by.y = "GeneMappingID", all.x = TRUE)
  
  names(Edge.source)[names(Edge.source)=="GeneSymbol"] <- "source.ID"
  
  Edge.target <- merge(Edge.source, NodeInfo[, c("GeneMappingID", "GeneSymbol")], by.x = "target", by.y = "GeneMappingID", all.x = TRUE)
  
  names(Edge.target)[names(Edge.target)=="GeneSymbol"] <- "target.ID"
  
  #Set column for pathway groupings
  NodeInfo$Group <- "Novel"
  
  NodeInfo$Group[NodeInfo$GeneSymbol %in% path1_pathway.genes.matrix] <- path1_name
  
  NodeInfo$Group[NodeInfo$GeneSymbol %in% path2_pathway.genes.matrix] <- path2_name        
  
  NodeInfo$Group[NodeInfo$GeneSymbol %in% path3_pathway.genes.matrix] <- path3_name
  
  #Merge Pathway, and TRIAGEhits to NodeInfo
  Scores_and_nodes <- merge(NodeInfo[, c("GeneMappingID", "GeneSymbol", "Group")],                      #Pairing up the "Node Info" with the gene info (such as gene symbol and groupings)
                            TRIAGEhits[, c("GeneSymbol", "EntrezID", "Pathway")], 
                            by.x = "GeneSymbol", by.y = "GeneSymbol", all.x = T)
  
  
  #Aggregate the edges to be summarised to each genemap ID                                            #pulling together all genes that interactect with a specifc gene and putting them in one row seprated by comma
  Edge_source_summary <- aggregate(target.ID ~ source.ID, data = Edge.target, paste, collapse = ", ")
  Edge_target_summary <- aggregate(source.ID ~ target.ID, data = Edge.target, paste, collapse = ", ")
  
  
  #Align column names and stack data frames                                                         #The analysis assumed directionality of interactions, but we're ignoring it here, so combining the "target" and Source" to one dataframe.
  colnames(Edge_target_summary) <- c("source.ID", "target.ID")
  Edge_summary_stacked <- rbind(Edge_source_summary, Edge_target_summary)
  
  #Aggregate stacked data to get unique values
  Edge_summary <- aggregate(target.ID ~ source.ID, data = Edge_summary_stacked, paste, collapse = ", ")
  
  #Update Names                                                                                     #Now a dataframe is being created where each gene selected by TRIAGE (or part of the highlighted groups) has a list "Ntwrk.all" that lists all other genes from TRAIGE that it is predicted to interact with.
  colnames(Edge_summary) <- c("GeneSymbol", "Ntwrk.all")
  
  #Merge with scores 
  Scores_nodes_and_edges <- merge(Scores_and_nodes, Edge_summary, by.x = "GeneSymbol", by.y = "GeneSymbol", all.x = T)
  #Create column with counts for interacting genes                                               #Here a counter is added, this counts for each gene how many genes (within the TRIAGE set) are predicted to have interactions with it, (this allows to then list your hits based on how many interactions it has with other hits)
  for (i in 1:length(Scores_nodes_and_edges$Ntwrk.all)) {
    Scores_nodes_and_edges$Allnet.count[i] <- ifelse(is.na(Scores_nodes_and_edges$Ntwrk.all[i]) == T, 0, 
                                                     length(unlist(strsplit(Scores_nodes_and_edges$Ntwrk.all[i], ", "))))
  }
  #Like the "counter" for all the network genes, now setting up seperate counter for each pahway of interest. Within pathways you want to seperate wther the gene it is interacting with is also a "hit" in TRIAGE or is it just a gene in that pathway that is not a hit.
  ########################################### Create Matrices of the genes within each group that are also hits in the screen
  #Pathway#1 Hit Genes 
  path1_hits.genes <- filter(Scores_nodes_and_edges, Group == path1_name, EntrezID %in% TRIAGEhits.matrix)
  
  path1_hits.genes.matrix <- matrix(path1_hits.genes$GeneSymbol)
  
  #Pathway#2 Hit Genes 
  
  path2_hits.genes <- filter(Scores_nodes_and_edges, Group == path2_name, EntrezID %in% TRIAGEhits.matrix)
  path2_hits.genes.matrix <- matrix(path2_hits.genes$GeneSymbol)
  
  #Pathway#3 Hit Genes 
  
  path3_hits.genes <- filter(Scores_nodes_and_edges, Group == path3_name, EntrezID %in% TRIAGEhits.matrix)
  path3_hits.genes.matrix <- matrix(path3_hits.genes$GeneSymbol)
  
  #Now for each "pathway" or group you do a seperate count, first pulling out the genes from the network column related to it, then counting how many of it it is.          
  ########################################### Create Lists and Rankings by group------path1_
  #path1_ Genes
  #Create Blank
  Scores_nodes_and_edges$Ntwrk.path1_ <- NA                                                         #First creating a blank column where to place these values
  group.list <- path1_pathway.genes.matrix                                                          #Defining the group of genes (in matrix form) from which the network list should be comprised of.
  
  #Extract genes from network interaction associted with genes from group                        #Here the formula looks for what genes in the "group" are in the list of genes of "Ntwrk.all" and formats it in a comma seperated way.
  for (i in 1:length(group.list)) {
    temp.name <- group.list[i]
    out <- grep(paste0("\\b", temp.name, "\\b"), Scores_nodes_and_edges[["Ntwrk.all"]], value = F, fixed = F)
    Scores_nodes_and_edges$Ntwrk.path1_[out] <- paste(temp.name, Scores_nodes_and_edges$Ntwrk.path1_[out], sep = ", ")
  }
  
  #remove "NA" from rows with genes in them
  Scores_nodes_and_edges$Ntwrk.path1_ <- gsub(", NA", "", Scores_nodes_and_edges$Ntwrk.path1_)
  
  #Create column with counts for interacting genes                                            #adding a count for how many interactions each gene has with gene in the "group"
  for (i in 1:length(Scores_nodes_and_edges$Ntwrk.path1_)) {
    Scores_nodes_and_edges$path1_net.count[i] <- ifelse(is.na(Scores_nodes_and_edges$Ntwrk.path1_[i]) == T, 0, 
                                                        length(unlist(strsplit(Scores_nodes_and_edges$Ntwrk.path1_[i], ", "))))
  }
  
  ########################################### Create Lists and Rankings by group----path2_          #Same process is repeated for genes in another group (there'e probably a less clunky way to do this, but I just copied over and over.)
  #path2_ Genes
  #Create Blank
  Scores_nodes_and_edges$Ntwrk.path2_ <- NA
  group.list <- path2_pathway.genes.matrix
  
  #Extract genes from network interaction associted with genes from group
  for (i in 1:length(group.list)) {
    temp.name <- group.list[i]
    out <- grep(paste0("\\b", temp.name, "\\b"), Scores_nodes_and_edges[["Ntwrk.all"]], value = F, fixed = F)
    Scores_nodes_and_edges$Ntwrk.path2_[out] <- paste(temp.name, Scores_nodes_and_edges$Ntwrk.path2_[out], sep = ", ")
  }
  
  #remove "NA" from rows with genes in them
  Scores_nodes_and_edges$Ntwrk.path2_ <- gsub(", NA", "", Scores_nodes_and_edges$Ntwrk.path2_)
  
  #Create column with counts for interacting genes
  for (i in 1:length(Scores_nodes_and_edges$Ntwrk.path2_)) {
    Scores_nodes_and_edges$path2_net.count[i] <- ifelse(is.na(Scores_nodes_and_edges$Ntwrk.path2_[i]) == T, 0, 
                                                        length(unlist(strsplit(Scores_nodes_and_edges$Ntwrk.path2_[i], ", "))))
  }
  
  ########################################### Create Lists and Rankings by group----path3_
  #path3_ Genes
  #Create Blank
  Scores_nodes_and_edges$Ntwrk.path3_ <- NA
  group.list <- path3_pathway.genes.matrix
  
  #Extract genes from network interaction associted with genes from group
  for (i in 1:length(group.list)) {
    temp.name <- group.list[i]
    out <- grep(paste0("\\b", temp.name, "\\b"), Scores_nodes_and_edges[["Ntwrk.all"]], value = F, fixed = F)
    Scores_nodes_and_edges$Ntwrk.path3_[out] <- paste(temp.name, Scores_nodes_and_edges$Ntwrk.path3_[out], sep = ", ")
  }
  
  #remove "NA" from rows with genes in them
  Scores_nodes_and_edges$Ntwrk.path3_ <- gsub(", NA", "", Scores_nodes_and_edges$Ntwrk.path3_)
  
  #Create column with counts for interacting genes
  for (i in 1:length(Scores_nodes_and_edges$Ntwrk.path3_)) {
    Scores_nodes_and_edges$path3_net.count[i] <- ifelse(is.na(Scores_nodes_and_edges$Ntwrk.path3_[i]) == T, 0, 
                                                        length(unlist(strsplit(Scores_nodes_and_edges$Ntwrk.path3_[i], ", "))))
  }
  
  ##################################### Create list and rankings by group hits ##############################             #This is the same process as above only instead of looking for the number of genes in the group it looks at genes in the group that are also "hits"
  
  
  ############################# path1_ Hits Genes
  
  #Create Blank
  Scores_nodes_and_edges$Ntwrk.path1_.hits <- NA
  group.list <- path1_hits.genes.matrix
  
  #Extract genes from network interaction associted with genes from group
  for (i in 1:length(group.list)) {
    temp.name <- group.list[i]
    out <- grep(paste0("\\b", temp.name, "\\b"), Scores_nodes_and_edges[["Ntwrk.all"]], value = F, fixed = F)
    Scores_nodes_and_edges$Ntwrk.path1_.hits[out] <- paste(temp.name, Scores_nodes_and_edges$Ntwrk.path1_.hits[out], sep = ", ")
  }
  
  #remove "NA" from rows with genes in them
  Scores_nodes_and_edges$Ntwrk.path1_.hits <- gsub(", NA", "", Scores_nodes_and_edges$Ntwrk.path1_)
  
  #Create column with counts for interacting genes
  for (i in 1:length(Scores_nodes_and_edges$Ntwrk.path1_.hits)) {
    Scores_nodes_and_edges$path1_hits.net.count[i] <- ifelse(is.na(Scores_nodes_and_edges$Ntwrk.path1_.hits[i]) == T, 0, 
                                                             length(unlist(strsplit(Scores_nodes_and_edges$Ntwrk.path1_[i], ", "))))
  }
  
  ############################# path2_ Hits Genes
  
  #Create Blank
  Scores_nodes_and_edges$Ntwrk.path2_.hits <- NA
  group.list <- path2_hits.genes.matrix
  
  #Extract genes from network interaction associted with genes from group
  for (i in 1:length(group.list)) {
    temp.name <- group.list[i]
    out <- grep(paste0("\\b", temp.name, "\\b"), Scores_nodes_and_edges[["Ntwrk.all"]], value = F, fixed = F)
    Scores_nodes_and_edges$Ntwrk.path2_.hits[out] <- paste(temp.name, Scores_nodes_and_edges$Ntwrk.path2_.hits[out], sep = ", ")
  }
  
  #remove "NA" from rows with genes in them
  Scores_nodes_and_edges$Ntwrk.path2_.hits <- gsub(", NA", "", Scores_nodes_and_edges$Ntwrk.path2_.hits)
  
  #Create column with counts for interacting genes
  for (i in 1:length(Scores_nodes_and_edges$Ntwrk.path2_.hits)) {
    Scores_nodes_and_edges$path2_hits.net.count[i] <- ifelse(is.na(Scores_nodes_and_edges$Ntwrk.path2_.hits[i]) == T, 0, 
                                                             length(unlist(strsplit(Scores_nodes_and_edges$Ntwrk.path2_.hits[i], ", "))))
  }
  
  ############################# path3_ Hits Genes
  
  #Create Blank
  Scores_nodes_and_edges$Ntwrk.path3_.hits <- NA
  group.list <- path3_hits.genes.matrix
  
  #Extract genes from network interaction associted with genes from group
  for (i in 1:length(group.list)) {
    temp.name <- group.list[i]
    out <- grep(paste0("\\b", temp.name, "\\b"), Scores_nodes_and_edges[["Ntwrk.all"]], value = F, fixed = F)
    Scores_nodes_and_edges$Ntwrk.path3_.hits[out] <- paste(temp.name, Scores_nodes_and_edges$Ntwrk.path3_.hits[out], sep = ", ")
  }
  
  #remove "NA" from rows with genes in them
  Scores_nodes_and_edges$Ntwrk.path3_.hits <- gsub(", NA", "", Scores_nodes_and_edges$Ntwrk.path3_.hits)
  
  #Create column with counts for interacting genes
  for (i in 1:length(Scores_nodes_and_edges$Ntwrk.path3_.hits)) {
    Scores_nodes_and_edges$path3_hits.net.count[i] <- ifelse(is.na(Scores_nodes_and_edges$Ntwrk.path3_.hits[i]) == T, 0, 
                                                             length(unlist(strsplit(Scores_nodes_and_edges$Ntwrk.path3_.hits[i], ", "))))
  }
  
  
  ############################## Network Hits Total                                                 #here a counter is added for how many "hits" that are in the selected groups it's predicted to interact with.
  Scores_nodes_and_edges <- Scores_nodes_and_edges %>%
    mutate(Total_Path_Hits.net.count = path1_hits.net.count + path3_hits.net.count + path2_hits.net.count)
  
  
  ###write file                                                                                   #Final File is created
  
  #Organize column order
  Scores_nodes_and_edges <- Scores_nodes_and_edges[, c("EntrezID", "GeneSymbol", "GeneMappingID"
                                                       , "Group", "Pathway"
                                                       , "Allnet.count", "Ntwrk.all"
                                                       ,"path1_hits.net.count" ,"Ntwrk.path1_.hits"
                                                       ,"path2_hits.net.count" ,"Ntwrk.path2_.hits"
                                                       ,"path3_hits.net.count" ,"Ntwrk.path3_.hits"
                                                       ,"Total_Path_Hits.net.count"
                                                       ,"path1_net.count", "Ntwrk.path1_"
                                                       ,"path2_net.count", "Ntwrk.path2_"
                                                       ,"path3_net.count", "Ntwrk.path3_"
  )]
  
  names(Scores_nodes_and_edges)[names(Scores_nodes_and_edges) == "Ntwrk.path1_"] <- paste0("Ntwrk.", path1_name)
  names(Scores_nodes_and_edges)[names(Scores_nodes_and_edges) == "path1_net.count"] <- paste0(path1_name, "_net.count")
  names(Scores_nodes_and_edges)[names(Scores_nodes_and_edges) == "Ntwrk.path2_"] <- paste0("Ntwrk.", path2_name)
  names(Scores_nodes_and_edges)[names(Scores_nodes_and_edges) == "path2_net.count"] <- paste0(path2_name, "_net.count")
  names(Scores_nodes_and_edges)[names(Scores_nodes_and_edges) == "Ntwrk.path3_"] <- paste0("Ntwrk.", path3_name)
  names(Scores_nodes_and_edges)[names(Scores_nodes_and_edges) == "path3_net.count"] <- paste0(path3_name, "_net.count")
  names(Scores_nodes_and_edges)[names(Scores_nodes_and_edges) == "Ntwrk.path1_.hits"] <- paste0("Ntwrk.", path1_name, ".hits")
  names(Scores_nodes_and_edges)[names(Scores_nodes_and_edges) == "path1_hits.net.count"] <- paste0(path1_name, ".hits.net.count")
  names(Scores_nodes_and_edges)[names(Scores_nodes_and_edges) == "Ntwrk.path2_.hits"] <- paste0("Ntwrk.", path2_name, ".hits")
  names(Scores_nodes_and_edges)[names(Scores_nodes_and_edges) == "path2_hits.net.count"] <- paste0(path2_name, ".hits.net.count")
  names(Scores_nodes_and_edges)[names(Scores_nodes_and_edges) == "Ntwrk.path3_.hits"] <- paste0("Ntwrk.", path3_name, ".hits")
  names(Scores_nodes_and_edges)[names(Scores_nodes_and_edges) == "path3_hits.net.count"] <- paste0("Ntwrk.", path3_name, ".hits.net.count")
  
  
  
  message(TRIAGE.output, "**")
  write.csv(Scores_nodes_and_edges, "Ranking_HumanTNFScreen.csv")
  
  ############################################################################### Add visualization ##############################################################################
  
  #############**************************************************################
  ###############################################################################
  #                 Hirachical Edge Bundling
  ###############################################################################
  
  
  
  #Set up dataframe for groupings (Loc)                                          #Assigning a number to each group (this will enable the grouping later)
  NodeInfo$Loc <- 4
  NodeInfo$Loc[NodeInfo$EntrezID %in% path3_nonhits.matrix] <-7
  NodeInfo$Loc[NodeInfo$EntrezID %in% path3_hits.matrix] <-3
  NodeInfo$Loc[NodeInfo$EntrezID %in% path2_nonhits.matrix] <-6
  NodeInfo$Loc[NodeInfo$EntrezID %in% path2_hits.matrix] <-2
  NodeInfo$Loc[NodeInfo$EntrezID %in% path1_nonhits.matrix] <-5
  NodeInfo$Loc[NodeInfo$EntrezID %in% path1_hits.matrix] <-1
  
  #Set up dataframe for IDs                                                      #now creating a name in the format of groupNumber.geneSymbol
  NodeInfo$ID <- paste(NodeInfo$Loc, NodeInfo$GeneSymbol, sep = ".")
  names(NodeInfo)[names(NodeInfo)== "GeneSymbol"] <- "key"
  
  #Move ID column first
  NodeInfo = NodeInfo[,c('ID', 'GeneMappingID', 'key', 'Loc')]
  
  #Set up rel file                                                              #The Hirarchical edge bundle package needs to dataframes, a NodeInfor with information about the nodes and a "rel" file about the relationships to be highilighted.
  rel.source <- merge(EdgeInfo, NodeInfo[, c("GeneMappingID", "ID", "Loc")], by.x = "source", by.y = "GeneMappingID", all.x = TRUE)      #To create the rel file the "EdgeInfo" file is combined with teh NodeInfo information
  names(rel.source)[names(rel.source)=="ID"] <- "source.ID"
  names(rel.source)[names(rel.source)=="Loc"] <- "Loc.source"
  
  rel.target <- merge(rel.source, NodeInfo[, c("GeneMappingID", "ID", "Loc")], by.x = "target", by.y = "GeneMappingID", all.x = TRUE)
  names(rel.target)[names(rel.target)=="ID"] <- "target.ID"
  names(rel.target)[names(rel.target)=="Loc"] <- "Loc.target"
  
  #Remove non hits                                                   #There are more connections than we want to visualize so here some connections are being removed to simplify
  rel.target.filter <- filter(rel.target, Loc.source !=  Loc.target) #Here connections within the same group are ignored (assuming that we are only interested in intergroup connections)
  rel.7 <- filter(rel.target, Loc.source != 4 & Loc.target != 4 &    #Here connections between Novel genes and other Novels genes are pulled out to be added back in later.
                    Loc.source != 3 & Loc.target != 3 &
                    Loc.source != 2 & Loc.target != 2 &
                    Loc.source != 1 & Loc.target != 1 )
  
  # Netgraph with 2nd dimension of connections
  rel2.7 <- filter(rel.target, Loc.source == 4 & Loc.target == 4)    #Here connections between Novel genes and other Novels genes are pulled out to be added back in later.
  rel.target <- rbind(rel.target.filter, rel.7)                      #Intra-connections of Novel genes are added to list of inter-group connections
  
  # For netgraph with 2nd dimension of connections
  rel2.target <- rbind(rel.target.filter, rel2.7)                      #Intra-connections of Novel genes are added to list of inter-group connections
  
  rel.target <- filter(rel.target,                                   #Connections to "non-hits" are removed
                       Loc.source != 7 &
                         Loc.target != 7 &
                         Loc.source != 5 &
                         Loc.target != 5 &
                         Loc.source != 6 &
                         Loc.target != 6 )
  
  # For netgraph with 2nd dimension of connections
  rel2.target <- filter(rel2.target,                                   #Connections to "non-hits" are removed
                        Loc.source != 7 &
                          Loc.target != 7 &
                          Loc.source != 5 &
                          Loc.target != 5 &
                          Loc.source != 6 &
                          Loc.target != 6 )
  #Filter Node info                                                 #"Non hits" are removed from the NodeInfo list as well.
  NodeInfo <- filter(NodeInfo, Loc != 7 & Loc != 5 & Loc != 6)
  
  # Fix directionality of connection for members of group seven    #Hirarchical edge bundling colors the interactions based on the source so here the columns are switched to get the result wanted (shoudl Novel-group interactions be the color of the group or the color of the "novel" gene category)
  rel.target.filter <- filter(rel.target, Loc.target != 4)
  # For netgraph with 2nd dimension of connections
  rel2.target.filter <- filter(rel2.target, Loc.target != 4)
  rel.7 <- filter(rel.target, Loc.target == 4)
  # For netgraph with 2nd dimension of connections
  rel2.7 <- filter(rel2.target, Loc.target == 4)
  
  #Switchnames for 1st dimension graph only
  names(rel.7)[names(rel.7)=="target.ID"] <- "target.ID2"
  names(rel.7)[names(rel.7)=="source.ID"] <- "source.ID2"
  names(rel.7)[names(rel.7)=="target"] <- "target2"
  names(rel.7)[names(rel.7)=="source"] <- "source2"
  names(rel.7)[names(rel.7)=="target.ID2"] <- "source.ID"
  names(rel.7)[names(rel.7)=="source.ID2"] <- "target.ID"
  names(rel.7)[names(rel.7)=="target2"] <- "source"
  names(rel.7)[names(rel.7)=="source2"] <- "target"
  
  #Switchnames for 2nd dimension graph
  names(rel2.7)[names(rel2.7)=="target.ID"] <- "target.ID2"
  names(rel2.7)[names(rel2.7)=="source.ID"] <- "source.ID2"
  names(rel2.7)[names(rel2.7)=="target"] <- "target2"
  names(rel2.7)[names(rel2.7)=="source"] <- "source2"
  names(rel2.7)[names(rel2.7)=="target.ID2"] <- "source.ID"
  names(rel2.7)[names(rel2.7)=="source.ID2"] <- "target.ID"
  names(rel2.7)[names(rel2.7)=="target2"] <- "source"
  names(rel2.7)[names(rel2.7)=="source2"] <- "target"
  
  rel.target <- rbind(rel.target.filter, rel.7)
  rel2.target <- rbind(rel2.target.filter, rel2.7)
  
  #Flip columns to get pathway colors as links
  names(rel.target)[names(rel.target)=="target.ID"] <- "target.ID2"
  names(rel.target)[names(rel.target)=="source.ID"] <- "source.ID2"
  names(rel.target)[names(rel.target)=="target"] <- "target2"
  names(rel.target)[names(rel.target)=="source"] <- "source2"
  names(rel.target)[names(rel.target)=="target.ID2"] <- "source.ID"
  names(rel.target)[names(rel.target)=="source.ID2"] <- "target.ID"
  names(rel.target)[names(rel.target)=="target2"] <- "source"
  names(rel.target)[names(rel.target)=="source2"] <- "target"
  
  rel <- rel.target[, c("source.ID", "target.ID")]
  names(rel)[names(rel)=="source.ID"] <- "V1"
  names(rel)[names(rel)=="target.ID"] <- "V2"
  
  #Flip columns to get pathway colors as links for 2nd dimension graph
  names(rel2.target)[names(rel2.target)=="target.ID"] <- "target.ID2"
  names(rel2.target)[names(rel2.target)=="source.ID"] <- "source.ID2"
  names(rel2.target)[names(rel2.target)=="target"] <- "target2"
  names(rel2.target)[names(rel2.target)=="source"] <- "source2"
  names(rel2.target)[names(rel2.target)=="target.ID2"] <- "source.ID"
  names(rel2.target)[names(rel2.target)=="source.ID2"] <- "target.ID"
  names(rel2.target)[names(rel2.target)=="target2"] <- "source"
  names(rel2.target)[names(rel2.target)=="source2"] <- "target"
  
  rel2 <- rel2.target[, c("source.ID", "target.ID")]
  names(rel2)[names(rel2)=="source.ID"] <- "V1"
  names(rel2)[names(rel2)=="target.ID"] <- "V2"
  
  # Remove nodes that do not have connection to the selected pathways
  rel.V1.matrix <- as.matrix((unique(rel$V1)))
  rel.V2.matrix <- as.matrix((unique(rel$V2)))
  rel.genes.matrix <- as.matrix(unique(rbind(rel.V1.matrix, rel.V2.matrix)))
  
  NodeInfo1 <<- filter(NodeInfo, Loc != 4 | ID %in% rel.genes.matrix)
  NodeInfo2 <<- NodeInfo
  
  #Generate the graph
  g <- graph.data.frame(rel, directed=T, vertices=NodeInfo1)
  g2 <- graph.data.frame(rel2, directed=T, vertices=NodeInfo2)
  
  clr <- as.factor(V(g)$Loc)
  clr2 <- as.factor(V(g2)$Loc)
  
  if(length(selectedRows) == 3){
    levels(clr) <- c("red", "darkblue", "saddlebrown", "green")  #Four colors are chosen since there are four groups including the other TRIAGE hit genes
    levels(clr2) <- c("red", "darkblue", "saddlebrown", "green")  #Four colors are chosen since there are four groups including the other TRIAGE hit genes
  }
  if(length(selectedRows) == 2){
    levels(clr) <- c("red", "darkblue", "green")  #Three colors are chosen since there are three groups
    levels(clr2) <- c("red", "darkblue", "green")  #Three colors are chosen since there are three groups
  }
  if(length(selectedRows) == 1){
    levels(clr) <- c("red", "green")  #Two colors are chosen since there are two groups
    levels(clr2) <- c("red", "green")  #Two colors are chosen since there are two groups
  }
  
  V(g)$color <- as.character(clr)
  V(g)$size = degree(g)*5
  
  V(g2)$color <- as.character(clr2)
  V(g2)$size = degree(g2)*5
  
  # igraph static plot
  #plot(g, layout = layout.circle, vertex.label=NA)
  
  Chimera1 <<- edgebundleR::edgebundle(g, tension = 0.8, fontsize = 8)       
  Chimera2 <<- edgebundleR::edgebundle(g2, tension = 0.8, fontsize = 3)       
  
  # Create 1st dimension networkD3 object
  g11 = g
  g11_wc <- cluster_walktrap(g11)
  g11_members <- membership(g11_wc)
  
  # Convert to object suitable for networkD3
  g11_d3 <<- igraph_to_networkD3(g11, group = g11_members)
  g11_vis <<- toVisNetworkData(g11)
    
  # Create 2nd dimension networkD3 object
  g22 = g2
  g22_wc <- cluster_walktrap(g22)
  g22_members <- membership(g22_wc)
  
  # Convert to object suitable for networkD3
  g22_d3 <<- igraph_to_networkD3(g22, group = g22_members)
  g22_vis <<- toVisNetworkData(g22)
  
  # Add a legend box on the html page
  if(length(selectedRows) == 3){
    graphLegend <<- sprintf('
  <div id="htmlwidget_container">
                            <form style="width: 360px; margin: 0 auto; color: grey;">
                            <fieldset>
                            <legend>Network Graph Colors:</legend>
                            <font color="red" face="courier"><b>&nbsp;&nbsp;Red:</b></font><font size="-1" color="red"> %s</font><br>
                            <font color="darkblue" face="courier"><b>&nbsp;Blue:</b></font><font size="-1" color="darkblue"> %s</font><br>
                            <font color="saddlebrown" face="courier"><b>Brown:</b></font><font size="-1"color="saddlebrown"> %s</font><br>
                            <font color="green" face="courier"><b>Green:</b></font><font size="-1" color="green"> %s</font><br>
                            </fieldset>
                            </form>',
                            path1_name, path2_name, path3_name, "other TRIAGE hit genes")
  }
  if(length(selectedRows) == 2){
    graphLegend <<- sprintf('
  <div id="htmlwidget_container">
                            <form style="width: 360px; margin: 0 auto; color: grey">
                            <fieldset>
                            <legend>Network Graph Colors:</legend>
                            <font color="red" face="courier"><b>&nbsp;&nbsp;Red:</b></font><font size="-1" color="red"> %s</font><br>
                            <font color="darkblue" face="courier"><b>&nbsp;Blue:</b></font><font size="-1" color="darkblue"> %s</font><br>
                            <font color="green" face="courier"><b>Green:</b></font><font size="-1" color="green"> %s</font><br>
                            </fieldset>
                            </form>',
                            path1_name, path2_name, "other TRIAGE hit genes")
  }
  if(length(selectedRows) == 1){
    graphLegend <<- sprintf('
  <div id="htmlwidget_container">
                            <form style="width: 360px; margin: 0 auto; color: grey;">
                            <fieldset>
                            <legend>Network Graph Colors:</legend>
                            <font color="red" face="courier"><b>&nbsp;&nbsp;Red:</b></font><font size="-1" color="red"> %s</font><br>
                            <font color="green" face="courier"><b>Green:</b></font><font size="-1" color="green"> %s</font><br>
                            </fieldset>
                            </form>',
                            path1_name, "other TRIAGE hit genes")
  }
  


  #Places (2) where plot will be saved to
  #setwd(TRIAGE.output)    
  #saveEdgebundle(Chimera1, "Chimera_STRINGHi_against_selectedPathways_1st.hits.html")
  #saveEdgebundle(Chimera2, "Chimera_STRINGHi_against_selectedPathways_2nd.hits.html")
  
  if(grepl('shiny', outputDir)){
    saveEdgebundle(Chimera1,file = "/srv/shiny-server/Chimera_STRINGHi_against_selectedPathways_1st.hits.html")
    saveEdgebundle(Chimera2,file = "/srv/shiny-server/Chimera_STRINGHi_against_selectedPathways_2nd.hits.html")
  }else{
    saveEdgebundle(Chimera1,file = "/Library/WebServer/Documents/Chimera_STRINGHi_against_selectedPathways_1st.hits.html")
    saveEdgebundle(Chimera2,file = "/Library/WebServer/Documents/Chimera_STRINGHi_against_selectedPathways_2nd.hits.html")
  }

  # # Add figure legend only if created when 1-3 pathways were selected
  # if(exists("graphLegend")){
  # 
  #   # Put a legend in the HTML file (inputOutput directory)
  #   #inHTML  <- readLines("Chimera_STRINGHi_MoTNF.hits.html")
  #   inHTML  <- readLines("Chimera_STRINGHi_against_selectedPathways_1st.hits.html")
  #   inHTML  <- readLines("Chimera_STRINGHi_against_selectedPathways_2nd.hits.html")
  #   outHTML  <- gsub(pattern = '<div id="htmlwidget_container">', replace = graphLegend, x = inHTML)
  #   writeLines(outHTML, con="Chimera_STRINGHi_against_selectedPathways_1st.hits.html")
  #   writeLines(outHTML, con="Chimera_STRINGHi_against_selectedPathways_2nd.hits.html")
  # 
  #   # Put a legend in the HTML file (localhost)
  #   if(grepl('shiny', outputDir)){
  #     inHTML2  <- readLines("/srv/shiny-server/Chimera_STRINGHi_against_selectedPathways_1st.hits.html")
  #   }else{
  #     inHTML2  <- readLines("/Library/WebServer/Documents/Chimera_STRINGHi_MoTNF.hits.html")
  #   }
  # 
  #   outHTML2  <- gsub(pattern = '<div id="htmlwidget_container">', replace = graphLegend, x = inHTML2)
  # 
  #   if(grepl('shiny', outputDir)){
  #     writeLines(outHTML2, con="/srv/shiny-server/Chimera_STRINGHi_against_selectedPathways_1st.hits.html")
  #     writeLines(outHTML2, con="/srv/shiny-server/Chimera_STRINGHi_against_selectedPathways_2nd.hits.html")
  #   }else{
  #     writeLines(outHTML2, con="/Library/WebServer/Documents/Chimera_STRINGHi_against_selectedPathways_1st.hits.html")
  #     writeLines(outHTML2, con="/Library/WebServer/Documents/Chimera_STRINGHi_against_selectedPathways_2nd.hits.html")
  #   }
  # }
  return(TRUE)
}
