# Create seperate files for all six cirtteria in STRING seperated by confidence level and species
# October 3, 2018
# Samuel Katz
########################

start.time <- Sys.time()

###Updated STRING Files, seperating direct interactin from interlogs (Ones based on orthologoues)
#Libraries 
library ('dplyr')
library('igraph')
library('data.table')

#Working Directories
DatabaseWD <- "~/Documents/TRIAGE/NetworkAnalysis/Update_Oct2018/CSVfiles/"
ConversionWD <- "~/Documents/IAM/STRINGnet/IDconversionFiles/"

#get .txt file
setwd(DatabaseWD)
String.human.full <- fread("9606.protein.links.full.v10.5.txt", header = T)
String.mouse.full <- fread("10090.protein.links.full.v10.5.txt", header = T)


# function to create full human and mouse dataframes
create_df <- function(df){
  splitIDs <- apply(df[, c('protein1', 'protein2')], 2, function(y) sub(".*\\.", "", y))
  df[, protein1 := splitIDs[, 1]]
  df[, protein2 := splitIDs[, 2]]
  setnames(df, c(1:2, 10), c('protID.1', 'protID.2', 'experimental'))
  return(df)
}

#Add columns seperating protein name from taxonomy code
##Human
String.human.full <- create_df(String.human.full)

##Mouse
String.mouse.full <- create_df(String.mouse.full)


##############Get EntrezIDs

#get conversion file
setwd(ConversionWD)
gene2ensembl <- fread("gene2ensembl_header.txt", header = T)

# function to create entrez full human and mouse dataframes
create_entrez <- function(df, stringDF, type='human'){
  code = ifelse(type == 'human', 9606, 10090)
  df_f = df[tax_id == code, c("GeneID", "Ensembl_protein_identifier")]
  m1 <- merge(stringDF, df_f, 
               by.x = "protID.1", by.y = "Ensembl_protein_identifier", 
               all.x = T, incomparables = 0)
  colnames(m1)[colnames(m1)=='GeneID'] = "entrez1"
  m2 <- merge(m1, df_f, 
              by.x = "protID.2", by.y = "Ensembl_protein_identifier", 
              all.x = T, incomparables = 0)
  colnames(m2)[colnames(m2)=='GeneID'] = "entrez2"
  df_entrez <- m2 %>% filter(!is.na(entrez1) & !is.na(entrez2))
  df_entrez = data.table(df_entrez)
  return(df_entrez)
}

## Human
String.human.full.entrez = create_entrez(gene2ensembl, String.human.full)

## Mouse
String.mouse.full.entrez = create_entrez(gene2ensembl, String.mouse.full, "mouse")

#######################
## Create igraphs

# function to assign subfiles based on confidence level
assign_conf <- function(j, crit, df){
  cols = c(crit, paste('entrez', c(1,2), sep=''))
  if(j == 'lowConf'){
    val = df[get(crit) >= 150 & get(crit) < 400, ..cols]
  } 
  else if(j == 'midConf'){
    val = df[get(crit) >= 400 & get(crit) < 700, ..cols]
  } 
  else{
    val = df[get(crit) >= 700, ..cols]
  }
  return(val)
}

# function to create a human/ mouse igraphs and save the data in wd
# create_igraph <- function(key, type, level){
#   types = c('human', 'mouse')
#   conf.levels = c('lowConf', 'midConf', 'highConf')
#   
#   for(type in types){
#     print(paste0(key, '-', type))
#     df.name = paste0('String.', type, '.full.entrez')
#     df = eval(parse(text=df.name))
#     
#     for(j in conf.levels){
#       val = assign_conf(j, key, df)
#       str.name = paste0("String.", type, ".", key, ".", j, ".igraph.Rdata")
#       i_graph = graph.data.frame(d = val, directed = FALSE)
#       return(i_graph)
#     }
#   }
# }

create_igraph <- function(key, type, level){
  print(paste0(key, '-', type, '-', level))
  df.name = paste0('String.', type, '.full.entrez')
  df = eval(parse(text=df.name))
  val = assign_conf(level, key, df)
  cols = paste('entrez', c(1,2), sep='')
  i_graph = graph.data.frame(d = val[, ..cols], directed = FALSE)
  E(i_graph)$weights = val[,get(key)]
  E(i_graph)$datasource = key
  return(i_graph)
}

create_special_igraphs <- function(str1, str2){
  load(paste0(str1, '.Rdata'))
  load(paste0(str2, '.Rdata'))
  G = graph.union(get(str1), get(str2))
  E(G)$weights = apply(cbind(E(G)$weights_1, E(G)$weights_2), 1, max, na.rm = T)
  datasources = apply(cbind(E(G)$weights_1, E(G)$weights_2), 1, which.max)
  E(G)$datasource = ifelse(datasources==1, E(G)$datasource_1, E(G)$datasource_2)
  return(G)
}

# to write into local file system
#setwd('/My Documents/Update_Oct2018/StringIgraphs/')

# to write into TRIAGE...
setwd("/Users/sakatz/TRIAGE_old2/app/data/Networks")

key = c('coexpression', 'cooccurence', 'database', 'experimental', 'fusion', 'neighborhood', 'textmining')

String.human.coexpression.highConf.igraph = create_igraph(key[1], 'human', 'highConf')
save(String.human.coexpression.highConf.igraph, file = 'String.human.coexpression.highConf.igraph.Rdata')
rm(String.human.coexpression.highConf.igraph)

String.human.coexpression.lowConf.igraph = create_igraph(key[1], 'human', 'lowConf')
save(String.human.coexpression.lowConf.igraph, file = 'String.human.coexpression.lowConf.igraph.Rdata')
rm(String.human.coexpression.lowConf.igraph)

String.human.coexpression.midConf.igraph = create_igraph(key[1], 'human', 'midConf')
save(String.human.coexpression.midConf.igraph, file = 'String.human.coexpression.midConf.igraph.Rdata')
rm(String.human.coexpression.midConf.igraph)

String.human.cooccurence.highConf.igraph = create_igraph(key[2], 'human', 'highConf')
save(String.human.cooccurence.highConf.igraph, file = 'String.human.cooccurence.highConf.igraph.Rdata')
rm(String.human.cooccurence.highConf.igraph)

String.human.cooccurence.lowConf.igraph = create_igraph(key[2], 'human', 'lowConf')
save(String.human.cooccurence.lowConf.igraph, file = 'String.human.cooccurence.lowConf.igraph')
rm(String.human.cooccurence.lowConf.igraph)

String.human.cooccurence.midConf.igraph = create_igraph(key[2], 'human', 'midConf')
save(String.human.cooccurence.midConf.igraph, file = 'String.human.cooccurence.midConf.igraph.Rdata')
rm(String.human.cooccurence.midConf.igraph)

String.human.database.highConf.igraph = create_igraph(key[3], 'human', 'highConf')
save(String.human.database.highConf.igraph, file = 'String.human.database.highConf.igraph.Rdata')
rm(String.human.database.highConf.igraph)

String.human.database.lowConf.igraph = create_igraph(key[3], 'human', 'lowConf')
save(String.human.database.lowConf.igraph, file = 'String.human.database.lowConf.igraph.Rdata')
rm(String.human.database.lowConf.igraph)

String.human.database.midConf.igraph = create_igraph(key[3], 'human', 'midConf')
save(String.human.database.midConf.igraph, file = 'String.human.database.midConf.igraph.Rdata')
rm(String.human.database.midConf.igraph)

String.human.experimental.highConf.igraph = create_igraph(key[4], 'human', 'highConf')
save(String.human.experimental.highConf.igraph, file = 'String.human.experimental.highConf.igraph.Rdata')
rm(String.human.experimental.highConf.igraph)

String.human.experimental.lowConf.igraph = create_igraph(key[4], 'human', 'lowConf')
save(String.human.experimental.lowConf.igraph, file = 'String.human.experimental.lowConf.igraph.Rdata')
rm(String.human.experimental.lowConf.igraph)

String.human.experimental.midConf.igraph = create_igraph(key[4], 'human', 'midConf')
save(String.human.experimental.midConf.igraph, file = 'String.human.experimental.midConf.igraph.Rdata')
rm(String.human.experimental.midConf.igraph)

String.human.fusion.highConf.igraph = create_igraph(key[5], 'human', 'highConf')
save(String.human.fusion.highConf.igraph, file = 'String.human.fusion.highConf.igraph.Rdata')
rm(String.human.fusion.highConf.igraph)

String.human.fusion.lowConf.igraph    = create_igraph(key[5], 'human', 'lowConf')
save(String.human.fusion.lowConf.igraph, file = 'String.human.fusion.lowConf.igraph.Rdata')
rm(String.human.fusion.lowConf.igraph)

String.human.fusion.midConf.igraph = create_igraph(key[5], 'human', 'midConf')
save(String.human.fusion.midConf.igraph, file = 'String.human.fusion.midConf.igraph.Rdata')
rm(String.human.fusion.midConf.igraph)

String.human.neighborhood.highConf.igraph = create_igraph(key[6], 'human', 'highConf')
save(String.human.neighborhood.highConf.igraph, file = 'String.human.neighborhood.highConf.igraph.Rdata')
rm(String.human.neighborhood.highConf.igraph)

String.human.neighborhood.lowConf.igraph = create_igraph(key[6], 'human', 'lowConf')
save(String.human.neighborhood.lowConf.igraph, file = 'String.human.neighborhood.lowConf.igraph.Rdata')
rm(String.human.neighborhood.lowConf.igraph)

String.human.neighborhood.midConf.igraph = create_igraph(key[6], 'human', 'midConf')
save(String.human.neighborhood.midConf.igraph, file = 'String.human.neighborhood.midConf.igraph.Rdata')
rm(String.human.neighborhood.midConf.igraph)

String.human.textmining.highConf.igraph = create_igraph(key[7], 'human', 'highConf')
save(String.human.textmining.highConf.igraph, file = 'String.human.textmining.highConf.igraph.Rdata')
rm(String.human.textmining.highConf.igraph)

String.human.textmining.lowConf.igraph = create_igraph(key[7], 'human', 'lowConf')
save(String.human.textmining.lowConf.igraph, file = 'String.human.textmining.lowConf.igraph.Rdata')
rm(String.human.textmining.lowConf.igraph)

String.human.textmining.midConf.igraph = create_igraph(key[7], 'human', 'midConf')
save(String.human.textmining.midConf.igraph, file = 'String.human.textmining.midConf.igraph.Rdata')
rm(String.human.textmining.midConf.igraph)

String.mouse.coexpression.highConf.igraph = create_igraph(key[1], 'mouse', 'highConf')
save(String.mouse.coexpression.highConf.igraph, file = 'String.mouse.coexpression.highConf.igraph.Rdata')
rm(String.mouse.coexpression.highConf.igraph)

String.mouse.coexpression.lowConf.igraph = create_igraph(key[1], 'mouse', 'lowConf')
save(String.mouse.coexpression.lowConf.igraph, file = 'String.mouse.coexpression.lowConf.igraph.Rdata')
rm(String.mouse.coexpression.lowConf.igraph)

String.mouse.coexpression.midConf.igraph = create_igraph(key[1], 'mouse', 'midConf')
save(String.mouse.coexpression.midConf.igraph, file = 'String.mouse.coexpression.midConf.igraph.Rdata')
rm(String.mouse.coexpression.midConf.igraph)

String.mouse.cooccurence.highConf.igraph = create_igraph(key[2], 'mouse', 'highConf')
save(String.mouse.cooccurence.highConf.igraph, file = 'String.mouse.cooccurence.highConf.igraph.Rdata')
rm(String.mouse.cooccurence.highConf.igraph)

String.mouse.cooccurence.lowConf.igraph =  create_igraph(key[2], 'mouse', 'lowConf')
save(String.mouse.cooccurence.lowConf.igraph, file = 'String.mouse.cooccurence.lowConf.igraph.Rdata')
rm(String.mouse.cooccurence.lowConf.igraph)

String.mouse.cooccurence.midConf.igraph = create_igraph(key[2], 'mouse', 'midConf')
save(String.mouse.cooccurence.midConf.igraph, file = 'String.mouse.cooccurence.midConf.igraph.Rdata')
rm(String.mouse.cooccurence.midConf.igraph)

String.mouse.database.highConf.igraph = create_igraph(key[3], 'mouse', 'highConf')
save(String.mouse.database.highConf.igraph, file = 'String.mouse.database.highConf.igraph.Rdata')
rm(String.mouse.database.highConf.igraph)

String.mouse.database.lowConf.igraph = create_igraph(key[3], 'mouse', 'lowConf')
save(String.mouse.database.lowConf.igraph, file = 'String.mouse.database.lowConf.igraph.Rdata')
rm(String.mouse.database.lowConf.igraph)

String.mouse.database.midConf.igraph = create_igraph(key[3], 'mouse', 'midConf')
save(String.mouse.database.midConf.igraph, file = 'String.mouse.database.midConf.igraph.Rdata')
rm(String.mouse.database.midConf.igraph)

String.mouse.experimental.highConf.igraph = create_igraph(key[4], 'mouse', 'highConf')
save(String.mouse.experimental.highConf.igraph, file = 'String.mouse.experimental.highConf.igraph.Rdata')
rm(String.mouse.experimental.highConf.igraph)

String.mouse.experimental.lowConf.igraph = create_igraph(key[4], 'mouse', 'lowConf')
save(String.mouse.experimental.lowConf.igraph, file = 'String.mouse.experimental.lowConf.igraph.Rdata')
rm(String.mouse.experimental.lowConf.igraph)

String.mouse.experimental.midConf.igraph = create_igraph(key[4], 'mouse', 'midConf')
save(String.mouse.experimental.midConf.igraph, file = 'String.mouse.experimental.midConf.igraph.Rdata')
rm(String.mouse.experimental.midConf.igraph)

String.mouse.fusion.highConf.igraph =      create_igraph(key[5], 'mouse', 'highConf')
save(String.mouse.fusion.highConf.igraph, file = 'String.mouse.fusion.highConf.igraph.Rdata')
rm(String.mouse.fusion.highConf.igraph)

String.mouse.fusion.lowConf.igraph = create_igraph(key[5], 'mouse', 'lowConf')
save(String.mouse.fusion.lowConf.igraph, file = 'String.mouse.fusion.lowConf.igraph.Rdata')
rm(String.mouse.fusion.lowConf.igraph)

String.mouse.fusion.midConf.igraph = create_igraph(key[5], 'mouse', 'midConf')
save(String.mouse.fusion.midConf.igraph, file = 'String.mouse.fusion.midConf.igraph.Rdata')
rm(String.mouse.fusion.midConf.igraph)

String.mouse.neighborhood.highConf.igraph = create_igraph(key[6], 'mouse', 'highConf')
save(String.mouse.neighborhood.highConf.igraph, file = 'String.mouse.neighborhood.highConf.igraph.Rdata')
rm(String.mouse.neighborhood.highConf.igraph)

String.mouse.neighborhood.lowConf.igraph = create_igraph(key[6], 'mouse', 'lowConf')
save(String.mouse.neighborhood.lowConf.igraph, file = 'String.mouse.neighborhood.lowConf.igraph.Rdata')
rm(String.mouse.neighborhood.lowConf.igraph)

String.mouse.neighborhood.midConf.igraph = create_igraph(key[6], 'mouse', 'midConf')
save(String.mouse.neighborhood.midConf.igraph, file = 'String.mouse.neighborhood.midConf.igraph.Rdata')
rm(String.mouse.neighborhood.midConf.igraph)

String.mouse.textmining.highConf.igraph =  create_igraph(key[7], 'mouse', 'highConf')
save(String.mouse.textmining.highConf.igraph, file = 'String.mouse.textmining.highConf.igraph.Rdata')
rm(String.mouse.textmining.highConf.igraph)

String.mouse.textmining.lowConf.igraph = create_igraph(key[7], 'mouse', 'lowConf')
save(String.mouse.textmining.lowConf.igraph, file = 'String.mouse.textmining.lowConf.igraph.Rdata')
rm(String.mouse.textmining.lowConf.igraph)

String.mouse.textmining.midConf.igraph = create_igraph(key[7], 'mouse', 'midConf')
save(String.mouse.textmining.midConf.igraph, file = 'String.mouse.textmining.midConf.igraph.Rdata')
rm(String.mouse.textmining.midConf.igraph)


# saving specialized files for quick retrieval of experimental and database

String.human.exp_and_data.highConf.igraph = create_special_igraphs("String.human.experimental.highConf.igraph", "String.human.database.highConf.igraph")
save(String.human.exp_and_data.highConf.igraph, file = 'String.human.exp_and_data.highConf.igraph.Rdata')
rm(String.human.exp_and_data.highConf.igraph)

String.human.exp_and_data.midConf.igraph = create_special_igraphs("String.human.experimental.midConf.igraph", "String.human.database.midConf.igraph")
save(String.human.exp_and_data.midConf.igraph, file = 'String.human.exp_and_data.midConf.igraph.Rdata')
rm(String.human.exp_and_data.midConf.igraph)

String.human.exp_and_data.lowConf.igraph = create_special_igraphs("String.human.experimental.lowConf.igraph", "String.human.database.lowConf.igraph")
save(String.human.exp_and_data.lowConf.igraph, file = 'String.human.exp_and_data.lowConf.igraph.Rdata')
rm(String.human.exp_and_data.lowConf.igraph)

String.mouse.exp_and_data.highConf.igraph = create_special_igraphs("String.mouse.experimental.highConf.igraph", "String.mouse.database.highConf.igraph")
save(String.mouse.exp_and_data.highConf.igraph, file = 'String.mouse.exp_and_data.highConf.igraph.Rdata')
rm(String.mouse.exp_and_data.highConf.igraph)

String.mouse.exp_and_data.midConf.igraph = create_special_igraphs("String.mouse.experimental.midConf.igraph", "String.mouse.database.midConf.igraph")
save(String.mouse.exp_and_data.midConf.igraph, file = 'String.mouse.exp_and_data.midConf.igraph.Rdata')
rm(String.mouse.exp_and_data.midConf.igraph)

String.mouse.exp_and_data.lowConf.igraph = create_special_igraphs("String.mouse.experimental.lowConf.igraph", "String.mouse.database.lowConf.igraph")
save(String.mouse.exp_and_data.lowConf.igraph, file = 'String.mouse.exp_and_data.lowConf.igraph.Rdata')
rm(String.mouse.exp_and_data.lowConf.igraph)



end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken












