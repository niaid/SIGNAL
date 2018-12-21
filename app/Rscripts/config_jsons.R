
# function to reclassify names as imports for two-way connection in network
change_name <- function(x){names(x) = "imports"; return(x)}

#function to configure dataframes for each json value
change_form <- function(l, i){data.frame(name = rep(names(l[i]), nrow(l[[i]])), l[[i]])}


config_json <- function(edges, dimension){
  
  # creates matrices of dimension name and edge matrix gene columns
  dim_m = cbind(dim1=dimension, edges[,1])
  dim_s = cbind(dimension, edges[,2])
  
  # collapses names and creates data frame of connections between level 1 and 2 genes
  name = apply(dim_m, 1, paste, collapse='.')
  imports = apply(dim_s, 1, paste, collapse='.')
  edge_df = data.frame(name, imports)
  
  # creates a larger list designating level 1 node connection to all other nodes
  group1 = lapply(split(edge_df, edge_df$name), `[`, 2)
  group2 = lapply(split(edge_df, edge_df$imports), `[`, 1)
  L = c(group1, group2)
  L = lapply(L, change_name)
  
  # json formation and output from list of dataframes
  df_L = list()
  for(i in 1:length(L)){df_L = append(df_L, list(change_form(L, i)))}
  json_1 = jsonlite::toJSON(df_L, 'columns')
  return(json_1)
}







