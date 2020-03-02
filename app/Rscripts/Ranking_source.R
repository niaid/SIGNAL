###################################################
## Setting up Network Analysis, jsons, and outputting
## TRIAGE hits
###################################################
# Kyle Webb, Sam Katz
# March 8, 2019
###################################################

Generate_NetworkGraph <- function(selectedRows, organism, G){
  #############################################################

  # TRIAGE hit data 
  TRIAGEhits <- TRIAGEoutput.condensed                                      

  #Get pathway file and pathway names (along with intersections)
  pathwayFile <- pathwayData    

  path1_name <<- sigPathways$Pathway[as.numeric(selectedRows[1])] 

  path1_pathway.genes <- filter(pathwayFile, PathwayName == path1_name)
  path1_pathway.genes.matrix <- matrix(path1_pathway.genes$EntrezID)
  
  path2_name <<- sigPathways$Pathway[as.numeric(selectedRows[2])]   
  path2_pathway.genes <- filter(pathwayFile, PathwayName == path2_name)
  path2_pathway.genes.matrix <- matrix(path2_pathway.genes$EntrezID)
  
  
  path3_name <<- sigPathways$Pathway[as.numeric(selectedRows[3])]    
  path3_pathway.genes <- filter(pathwayFile, PathwayName == path3_name)
  path3_pathway.genes.matrix <- matrix(path3_pathway.genes$EntrezID)

  
  path12_name = paste0(path1_name, ' & ', path2_name)
  path13_name = paste0(path1_name, ' & ', path3_name)
  path23_name = paste0(path2_name, ' & ', path3_name)
  path123_name = paste0(path1_name, ' & ', path2_name, ' & ', path3_name)

  # filtering pathway hits which are in the TRIAGE hits data
  path1_hits <- filter(TRIAGEhits, EntrezID %in% path1_pathway.genes.matrix)
  path1_hits.matrix <- matrix(path1_hits$EntrezID)
  
  path2_hits <- filter(TRIAGEhits, EntrezID %in% path2_pathway.genes.matrix)
  path2_hits.matrix <- matrix(path2_hits$EntrezID)
  
  path3_hits <- filter(TRIAGEhits, EntrezID %in% path3_pathway.genes.matrix)
  path3_hits.matrix <- matrix(path3_hits$EntrezID)
  
  #######################################                                            
  # redoing the pathway gene matrices to be of the gene symbol instead of the Entrez ID
  #Get Matrix of genes for each pathway in KEGG
  path1_pathway.genes.matrix <- matrix(path1_pathway.genes$GeneSymbol)
  
  path2_pathway.genes.matrix <- matrix(path2_pathway.genes$GeneSymbol)
  
  path3_pathway.genes.matrix <- matrix(path3_pathway.genes$GeneSymbol)
  
  TRIAGEhits.matrix <- matrix(TRIAGEhits$EntrezID)                      
  # Created a matrix of all the genes that are "hits"
  cat(file=stderr(), "Inside Generate_NetworkGraph64!\n")
  
  #Add GeneSymbols to Edge dataframe                                              
  Edge.source <- merge(EdgeInfo, NodeInfo[, c("GeneMappingID", "GeneSymbol")], by.x = "source", by.y = "GeneMappingID", all.x = TRUE)
  
  names(Edge.source)[names(Edge.source)=="GeneSymbol"] <- "source.ID"
  
  Edge.target <- merge(Edge.source, NodeInfo[, c("GeneMappingID", "GeneSymbol")], by.x = "target", by.y = "GeneMappingID", all.x = TRUE)
  
  names(Edge.target)[names(Edge.target)=="GeneSymbol"] <- "target.ID"

  
  #Set column for pathway groupings and colors for pathways
  colorMap = c("blue", #Pathway 1
               "red", #Pathway 2
               "saddlebrown", #Pathway 3
               "#a09d01", # 1 & 2
               "#ce6702", # 2 & 3
               "#6b5b3e", # 1 & 3
               "#6d1c8e", # 1 & 2 & 3
               "green" #Novel
  )
  
  NodeInfo$Group <- "Novel"
  NodeInfo$Color = last(colorMap)
  
  NodeInfo$Group[NodeInfo$GeneSymbol %in% path1_pathway.genes.matrix] <- path1_name
  NodeInfo$Color[NodeInfo$GeneSymbol %in% path1_pathway.genes.matrix] = colorMap[1]
  
  NodeInfo$Group[NodeInfo$GeneSymbol %in% path2_pathway.genes.matrix] <- path2_name 
  NodeInfo$Color[NodeInfo$GeneSymbol %in% path2_pathway.genes.matrix] = colorMap[2]
  
  NodeInfo$Group[NodeInfo$GeneSymbol %in% path3_pathway.genes.matrix] <- path3_name
  NodeInfo$Color[NodeInfo$GeneSymbol %in% path3_pathway.genes.matrix] = colorMap[3]
  
  # checks for hits available in intersection pathways
  path12_check = NodeInfo$GeneSymbol %in% path1_pathway.genes.matrix &
                  NodeInfo$GeneSymbol %in% path2_pathway.genes.matrix
  path13_check = NodeInfo$GeneSymbol %in% path1_pathway.genes.matrix &
                  NodeInfo$GeneSymbol %in% path3_pathway.genes.matrix
  path23_check = NodeInfo$GeneSymbol %in% path2_pathway.genes.matrix &
                  NodeInfo$GeneSymbol %in% path3_pathway.genes.matrix
  path123_check = NodeInfo$GeneSymbol %in% path1_pathway.genes.matrix &
                  NodeInfo$GeneSymbol %in% path2_pathway.genes.matrix &
                  NodeInfo$GeneSymbol %in% path3_pathway.genes.matrix
  
  # empty vector to be populated with hits that intersect multiple pathways
  extra_path_names = c()
  
  if(any(path12_check)){
    NodeInfo$Group[path12_check] <- path12_name
    NodeInfo$Color[path12_check] = colorMap[4]
    path12_name = paste0(path1_name, ' & ', path2_name)
    extra_path_names = append(extra_path_names, path12_name)
  }
  if(any(path23_check)){
    NodeInfo$Group[path23_check] <- path23_name
    NodeInfo$Color[path23_check] = colorMap[5]
    path23_name = paste0(path2_name, ' & ', path3_name)
    extra_path_names = append(extra_path_names, path23_name)
  }
  if(any(path13_check)){
    NodeInfo$Group[path13_check] <- path13_name
    NodeInfo$Color[path13_check] = colorMap[6]
    path13_name = paste0(path1_name, ' & ', path3_name)
    extra_path_names = append(extra_path_names, path13_name)
  }
  if(any(path123_check)){
    NodeInfo$Group[path123_check] <- path123_name
    NodeInfo$Color[path123_check] = colorMap[7]
    path123_name = paste0(path1_name, ' & ', path2_name, ' & ', path3_name)
    extra_path_names = append(extra_path_names, path123_name)
  }
  
  #Merge Pathway, and TRIAGEhits to NodeInfo
  Scores_and_nodes <- merge(NodeInfo[, c("GeneMappingID", "GeneSymbol", "Group")],                      
                            #Pairing up the "Node Info" with the gene info (such as gene symbol and groupings)
                            TRIAGEhits[, c("GeneSymbol", "EntrezID", "InputCategory", "Pathway")], 
                            by.x = "GeneSymbol", by.y = "GeneSymbol", all.x = T)
  
  
  #Aggregate the edges to be summarised to each genemap ID                                            
  #pulling together all genes that interactect with a specifc gene and putting them in one row seprated by comma
  Edge_source_summary <- aggregate(target.ID ~ source.ID, data = Edge.target, paste, collapse = ", ")
  Edge_target_summary <- aggregate(source.ID ~ target.ID, data = Edge.target, paste, collapse = ", ")
  
  
  #Align column names and stack data frames                                                         
  #The analysis assumed directionality of interactions, but we're ignoring it here, so combining the "target" and Source" to one dataframe.
  colnames(Edge_target_summary) <- c("source.ID", "target.ID")
  Edge_summary_stacked <- rbind(Edge_source_summary, Edge_target_summary)
  
  #Aggregate stacked data to get unique values
  Edge_summary <- aggregate(target.ID ~ source.ID, data = Edge_summary_stacked, paste, collapse = ", ")
  
  #Update Names                                                                                     
  #Now a dataframe is being created where each gene selected by TRIAGE (or part of the highlighted groups) has a list "Ntwrk.all" that lists all other genes from TRAIGE that it is predicted to interact with.
  colnames(Edge_summary) <- c("GeneSymbol", "Ntwrk.all")
  
  #Merge with scores 
  Scores_nodes_and_edges <- merge(Scores_and_nodes, Edge_summary, by.x = "GeneSymbol", by.y = "GeneSymbol", all.x = T)
  
  # Create column in Scores_nodes_and_edges with counts for intersecting pathways                                               #Here a counter is added, this counts for each gene how many genes (within the TRIAGE set) are predicted to have interactions with it, (this allows to then list your hits based on how many interactions it has with other hits)
  for (i in 1:length(Scores_nodes_and_edges$Ntwrk.all)) {
    Scores_nodes_and_edges$Allnet.count[i] <- ifelse(is.na(Scores_nodes_and_edges$Ntwrk.all[i]) == T, 0, 
                                                     length(unlist(strsplit(Scores_nodes_and_edges$Ntwrk.all[i], ", "))))
  }
  
  # organizing column order
  # Scores_nodes_and_edges <- Scores_nodes_and_edges[, c("EntrezID", "GeneSymbol", "GeneMappingID", "ConfidenceCategory"
  #                                                      , "Group", "Pathway"
  #                                                      , "Allnet.count", "Ntwrk.all")]
  
  Scores_nodes_and_edges <- Scores_nodes_and_edges[, c("EntrezID", "GeneSymbol", "InputCategory"
                                                       , "Group", "Pathway"
                                                       , "Allnet.count", "Ntwrk.all")]
  
  # renaming some of these columns...
  Scores_nodes_and_edges = Scores_nodes_and_edges %>% 
    rename("Pathway" = "TRIAGEpathways",
           "Allnet.count" = "Total.Network.Hits",
           "Ntwrk.all" = "All.Network.Hits")
  
  # creating vector for all pathway names for gene hits in multiple pathways
  
  if(length(extra_path_names)){
    all_pathway_names <<- c(sigPathways$Pathway[as.numeric(selectedRows)], extra_path_names)  
  }  else{
    all_pathway_names <<- sigPathways$Pathway[as.numeric(selectedRows)]
  }
  
  # Function and variables used in create_nodes_and_edges below for writing csv file name
  fix.names <- function(x){
    f.names = c()
    for(i in x){
      f.names = append(f.names, 
                       gsub("[[:space:] ]", "_", gsub("[^[:alnum:] ]", "", i)))
    }
    paste0(f.names, collapse="_")
  }
  effs = sigPathways$Pathway[as.numeric(selectedRows)]
  t.file.name = fix.names(effs)
  
  # Function to create final form of Scores_nodes_and_edges

  create_nodes_and_edges <- function(pn){
    n = length(pn)
    for(j in 1:n){
      
      i = pn[j]
      hit.genes <- filter(Scores_nodes_and_edges, Group == i & EntrezID %in% TRIAGEhits.matrix)
      hit.genes.matrix <- matrix(hit.genes$GeneSymbol)  
      temp.hits = 'hits'
      temp.counts = 'counts'
      Scores_nodes_and_edges[, temp.hits] <- NA
      Scores_nodes_and_edges[, temp.counts] <- 0
      group.list <- hit.genes.matrix
      all.out = c()
      
      for (k in 1:length(group.list)) {
        temp.name <- group.list[k]
        out <- grep(paste0("\\b", temp.name, "\\b"), Scores_nodes_and_edges[["All.Network.Hits"]], value = F, fixed = F)
        all.out = append(all.out, out)
        Scores_nodes_and_edges[out, temp.hits] <- paste(temp.name, Scores_nodes_and_edges[out, temp.hits], sep = ", ")
      }
      
      #remove "NA" from rows with genes in them
      Scores_nodes_and_edges[,temp.hits] <- gsub(", NA", "", Scores_nodes_and_edges[, temp.hits])
      
      for(q in all.out){
        Scores_nodes_and_edges[q, temp.counts] = length(unlist(strsplit(Scores_nodes_and_edges[q, temp.hits], ", ")))
      }
      
      # flipping the position of pathway hits and counts in the dataframe
      Scores_nodes_and_edges[, c(temp.hits, temp.counts)] = Scores_nodes_and_edges[, c(temp.counts, temp.hits)]
      # renaming columns for pathway names
      N = ncol(Scores_nodes_and_edges)
      sne.col.hits = paste0('Network.', i, '.Hits')
      sne.col.counts = paste0('Total.Network.', i)
      colnames(Scores_nodes_and_edges)[c(N-1, N)] = c(sne.col.counts, sne.col.hits)
      
    }
    
    # Group names that are "Novel" will be renamed to "Additional TRIAGE hits"
    Scores_nodes_and_edges = Scores_nodes_and_edges %>% 
      mutate(Group = ifelse(Group == 'Novel', 'Additional TRIAGE hits', Group))
    
    #here a counter is added for how many "hits" that are in the selected groups it's shown to interact with.
    sumcols = grep("Total.Network.", colnames(Scores_nodes_and_edges))
    
    if(n != 1){
      Scores_nodes_and_edges$Total.Network.Selected.Pathways = rowSums(Scores_nodes_and_edges[,sumcols])  
    }
    
    ### write file
    file.name = paste0(inputFilePrefix, "_TRIAGEnetwork_", t.file.name, ".csv")
    
    ret.list = list(Scores_nodes_and_edges, file.name)
    
    return(ret.list)
  }
  
  SNE_output = create_nodes_and_edges(all_pathway_names)
  
  Scores_nodes_and_edges <<- SNE_output[[1]]
  
  # scores nodes and edges filename
  file.name.snae <<- SNE_output[[2]]
  
  fwrite(Scores_nodes_and_edges, file.name.snae)
  
  ############################################################################### Add visualization ##############################################################################
  
  #############**************************************************################
  ###############################################################################
  #       Creating Node and Edge info for 1st and 2nd degree Networks
  ###############################################################################
  
  #Set up dataframe for IDs                                                      
  #now creating a name in the format of groupNumber.geneSymbol
  names(NodeInfo)[names(NodeInfo)== "GeneSymbol"] <- "key"
  names(NodeInfo)[names(NodeInfo)== "InputCategory"] = 'Confidence'
  
  NodeInfo$ID <- paste(NodeInfo$Group, NodeInfo$key, sep = ".")
  
  #Move ID column first
  NodeInfo = NodeInfo[,c('ID', 'GeneMappingID', 'key', 'Group', 'Confidence', 'Pathway', 'Color')]
  
  #Set up rel file                                                              
  #The Hirarchical edge bundle package needs to dataframes, a NodeInfor with information about the nodes and a "rel" file about the relationships to be highilighted.
  rel.source <- merge(EdgeInfo, NodeInfo[, c("GeneMappingID", "ID", "Group")], by.x = "source", by.y = "GeneMappingID", all.x = TRUE)      #To create the rel file the "EdgeInfo" file is combined with teh NodeInfo information
  names(rel.source)[names(rel.source)=="ID"] <- "source.ID"
  names(rel.source)[names(rel.source)=="Group"] <- "Group.source"
  
  rel.target <- merge(rel.source, NodeInfo[, c("GeneMappingID", "ID", "Group")], by.x = "target", by.y = "GeneMappingID", all.x = TRUE)
  names(rel.target)[names(rel.target)=="ID"] <- "target.ID"
  names(rel.target)[names(rel.target)=="Group"] <- "Group.target"
  
  # Remove hits within the same pathway (assuming only interested in inter-pathway interaction)                                  
  rel1.target <- filter(rel.target, Group.source !=  Group.target) 
  
  # Filter funciton that only includes 2nd degree novel genes IF
  #  they are connected to at least one novel gene that maps back to a pathway
  #  (i.e. Novel -> Novel -> Pathway)
  degree2.filter <- function(rel){
    markers = c()
    new.rel = filter(rel, Group.source == "Novel" & Group.target == "Novel")
    for(i in 1:nrow(new.rel)){
      r = new.rel[i,]
      st.filter = filter(rel, source.ID == r$source.ID | target.ID == r$target.ID)
      if(any(st.filter$Group.target!="Novel") || any(st.filter$Group.source!="Novel")){
        markers = append(markers, i)
      }
    }
    return(new.rel[markers,])
  }
  
  rel2.target = rbind(rel1.target, degree2.filter(rel.target))
  
  # labeling columns and adding weights and datasources for each network edge
  rel1 <- rel1.target[, c("source.ID", "target.ID")]
  names(rel1)[names(rel1)=="source.ID"] <- "V1"
  names(rel1)[names(rel1)=="target.ID"] <- "V2"
  rel1$weights = rel1.target$weights
  rel1$datasource = rel1.target$datasource
  
  rel2 <- rel2.target[, c("source.ID", "target.ID")]
  names(rel2)[names(rel2)=="source.ID"] <- "V1"
  names(rel2)[names(rel2)=="target.ID"] <- "V2"
  rel2$weights = rel2.target$weights
  rel2$datasource = rel2.target$datasource
  
  # Remove nodes that do not have connection to the selected pathways
  rel1.V1.matrix <- as.matrix((unique(rel1$V1)))
  rel1.V2.matrix <- as.matrix((unique(rel1$V2)))
  rel1.genes.matrix <- as.matrix(unique(rbind(rel1.V1.matrix, rel1.V2.matrix)))
 
  NodeInfo2 <<- NodeInfo
  NodeInfo1 <<- filter(NodeInfo, Group != "Novel" | ID %in% rel1.genes.matrix)
  
  #Generate the igraphs
  g <<- graph.data.frame(rel1, directed=T, vertices=NodeInfo1)
  g2 <<- graph.data.frame(rel2, directed=T, vertices=NodeInfo2)
  
  # Error message to make sure the generated graph is full
  if(length(E(g))==0 | length(V(g))==0){
    showModal(modalDialog(title="Warning:", HTML("<h3><font color=red>Criteria produced empty network. Session will restart.</font><h3>"),
                          easyClose = FALSE))
    Sys.sleep(5)
    session$reload()
  }
  
  
  ###############################################################################
  #                 Creating JSONs for D3
  ###############################################################################
  
  
  g11_vis <<- toVisNetworkData(g)

  # sources functions from config_jsons.R to build list (json_1df)
  json_1df <<- config_df(g11_vis$nodes, g11_vis$edges, all_pathway_names)
  json_1 = jsonlite::toJSON(json_1df, 'columns')

  # sends json_1 to custom_network.js file for first degree visualization
  session$sendCustomMessage(type="jsondata1",json_1)
  
  g22_vis <<- toVisNetworkData(g2)
  
  # sources functions from config_jsons.R to build list (json_2df)
  json_2df <<- config_df(g22_vis$nodes, g22_vis$edges, all_pathway_names)
  json_2 <- jsonlite::toJSON(json_2df, 'columns')
  
  # sends json_2 to custom_network2.js file for second degree visualization
  session$sendCustomMessage(type="jsondata2",json_2)
  
  # sources functions from config_jsons.R to build list for only Novel connections
  # nodes_3 = filter(g22_vis$nodes, Group=="Novel")
  # edges_3 = g22_vis$edges
  # edges_3 = filter(edges_3, gsub("\\..*","",edges_3$from) == "Novel" & 
  #                    gsub("\\..*","",edges_3$to) == "Novel")
  # json_3df <<- config_df(nodes_3, edges_3, all_pathway_names, onlyNovel=T)
  # json_3 <- jsonlite::toJSON(json_3df, 'columns')
  
  # sends json_2 to custom_network2.js file for second degree visualization
  # session$sendCustomMessage(type="jsondata3",json_3)
  
  PathNetName.output <<- paste0("PathNet_", inputFilePrefix, "_", t.file.name)
  
  return(TRUE)
}
