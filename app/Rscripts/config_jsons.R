
# function to reclassify names as imports for two-way connection in network
#change_name <- function(x){names(x) = "imports"; return(x)}

#function to configure dataframes for each json value
change_form <- function(l, i, nodes){
  key.name = gsub(".*\\.","",names(l[i]))
  node.info = nodes[which(nodes$key == key.name),]

  data.frame('name' = rep(names(l[i]), nrow(l[[i]])),
             'imports' = l[[i]][,1],
             'weights' = l[[i]][,2],
             'datasource' = l[[i]][,3],
             'Confidence' = rep(node.info$Confidence, nrow(l[[i]])))
             #'Pathway' = rep(node.info$Pathway, nrow(l[[i]])))
}

# function to combine names of edge dataframe
join_str <- function(str1, str2, names){
  c1 = substr(str1, 1, 1)
  c2 = substr(str2, 1, 1)
  val1 = as.numeric(c1)
  val2 = as.numeric(c2)
  st1 = substr(str1, 2, nchar(str1))
  st2 = substr(str2, 2, nchar(str2))
  if(val2 == 4){
    if(val1 != 4){
      c(paste0(names[val1], st1), paste0('Novel Genes', st2))  
    }
    else{
      c(paste0('Novel Genes', st1), paste0('Novel Genes', st2))
    }
  }
  else{
    c(paste0(names[val1], st1), paste0(names[val2], st2))
  }
}

#function to assign correct path name variables
assign_names <- function(names, edges){
  for(i in 1:nrow(edges)){
    edges[i,1:2] = join_str(edges[i,1], edges[i,2], names)
  }
  return(edges)
}

#function to merge the two lists of names and imports by the same name
merge_lists <- function(l.1, l.2){
  keys = unique(c(names(l.1), names(l.2)))
  L = lapply(keys, function(key){
    #print(key)
    l.1.col = l.1[key][[1]]
    l.2.col = l.2[key][[1]]
    check = any(is.null(l.1.col), is.null(l.2.col))
    if(!check){
      rbind(l.1.col, l.2.col)
    }
    else{
      w.full = ifelse(!is.null(l.1.col), 1, 2)
      list(l.1.col, l.2.col)[[w.full]]
    }
  })
  names(L) = keys
  return(L)
}

conf.f = function(x){
  x = substr(x, 1, nchar(x)-4)
  x = ifelse(x=='Med', 'Medium', x)
  return(x)
}

config_df <- function(nodes, edges, dimNames){
  #fix values of confidence levels
  nodes$Confidence = conf.f(nodes$Confidence)
  
  #assigns correct names of network pathways selected
  edges = assign_names(dimNames, edges)
  
  # change names of from and to to name and imports
  colnames(edges)[1:2] = c('name', 'imports')
  
  # creates a larger list designating level 1 node connection to all other nodes
  group1 = lapply(split(edges, edges$name), `[`, 2:ncol(edges))
  group2 = lapply(split(edges, edges$imports), `[`, c(1,3:ncol(edges)))
  group2 = lapply(group2, setNames, c('imports', 'weights', 'datasource')) #, 'keggConf', 'netConf', 'Pathway'))
  L = merge_lists(group1, group2)
  
  # json formation and output from list of dataframes
  df_L = list()
  for(i in 1:length(L)){df_L = append(df_L, list(change_form(L, i, nodes)))}
  return(df_L)
}

config_json <- function(nodes, edges, dimNames){
  df_L = config_df(nodes, edges, dimNames)
  json = jsonlite::toJSON(df_L, 'columns')
  
  return(json)
}

config_json2 <- function(nodes1, edges1, nodes2, edges2, dimNames){
  df1 = config_df(nodes1, edges1, dimNames)
  df2 = config_df(nodes2, edges2, dimNames)
  df_full = list(df1, df2)
  #df_full = df_full[-which(duplicated(df_full))]
  json = jsonlite::toJSON(df_full, 'columns')
  
  return(json)
}





