
# function to configure dataframes from merged lists of node-edge interaction
change_form <- function(l, i, nodes){
  key.name = gsub(".*\\.","",names(l[i]))
  node.info = nodes[which(nodes$key == key.name),]

  data.frame('name' = rep(names(l[i]), nrow(l[[i]])),
             'imports' = l[[i]][,1],
             'weights' = l[[i]][,2],
             'datasource' = l[[i]][,3],
             'Confidence' = rep(node.info$Confidence, nrow(l[[i]])),
             'order' = rep(i, nrow(l[[i]])),
             'Color' = rep(node.info$Color, nrow(l[[i]])))
             # 'Pathway' = rep(node.info$Pathway, nrow(l[[i]]))
}

# function to merge the two lists of names and imports by the same pathway.gene name
merge_lists <- function(l.1, l.2){
  keys = unique(c(names(l.1), names(l.2)))
  L = lapply(keys, function(key){
    l.1.col = l.1[key][[1]]
    # if(!is.null(l.1.col)) l.1.col = arrange(l.1.col, imports)
    l.2.col = l.2[key][[1]]
    # if(!is.null(l.2.col)) l.2.col = arrange(l.2.col, imports)
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

# Fixes "Medium" text for confidence values
conf.f <- function(x){
  x = substr(x, 1, nchar(x)-4)
  x = ifelse(x=='Med', 'Medium', x)
  return(x)
}

# Sorts the unique pathway names in the order of 
# Pathway 13, Pathway 1, Pathway 12, Pathway 2, Pathway 23, Pathway 3, Pathway 123, Novel Genes
orderedL.names <- function(uni.names){
  check = sapply(uni.names, function(s){grepl(' & ', s)})
  noAnds = which(sapply(uni.names, function(s){!grepl(' & ', s)}))
  double.check = sapply(uni.names, function(s){
    ifelse(str_count(s, " & ")==2, T, F)
  })
  n.uni.names = uni.names[1]
  if(any(check)){
    if(max(noAnds) == 2){
      n.uni.names = append(n.uni.names, c(names(which(check)), uni.names[2]))
    } 
    else{
      w.names = names(which(check))
      c.12 = paste0(uni.names[1], ' & ', uni.names[2]) %in% uni.names
      c.23 = paste0(uni.names[2], ' & ', uni.names[3]) %in% uni.names
      c.13 = paste0(uni.names[1], ' & ', uni.names[3]) %in% uni.names
      if(c.12){
        n.uni.names = append(n.uni.names, c(w.names[1], uni.names[2]))
      }
      else{
        n.uni.names = append(n.uni.names, uni.names[2])
      }
      if(c.23){
        w.23 = w.names[which(w.names == paste0(uni.names[2], ' & ', uni.names[3]))]
        n.uni.names = append(n.uni.names, c(w.23, uni.names[3]))
      }
      else{
        n.uni.names = append(n.uni.names, uni.names[3])
      }
      if(c.13){
        w.13 = w.names[which(w.names == paste0(uni.names[1], ' & ', uni.names[3]))]
        n.uni.names = append(n.uni.names, w.13)
      }
    }
  }
  else{
    n.uni.names = uni.names
  }
  if(any(double.check)){
    dd.c = which(double.check)
    n.uni.names = append(n.uni.names, names(double.check)[dd.c], after=0)
  }
  n.uni.names = append(n.uni.names, "Novel")
  return(n.uni.names)
}
  
# attributes the unique ordered pathway names to the merged list of node-edge interaction  
orderL <- function(L, dn){
  L.names = gsub("\\..*", "", names(L))
  check = sapply(dn, function(s){grepl(' & ', s)})
  if(!any(check)){
    ref = c(dn, "Novel")
  }
  else{
    ref = orderedL.names(dn)
  }
  if(length(ref) > 2){
    med.ref = ceiling(median(1:(length(ref)-1)))
    ref = ref[c(med.ref:length(ref), 1: (med.ref-1))]
  }
  ord.vect = as.vector(unlist(sapply(ref, function(s){
      brokenUp = which(L.names == s)
      match(sort(names(L[brokenUp])), names(L))
  })))
  l = L[ord.vect]
  names(l) = names(L)[ord.vect]
  return(l)
}


config_df <- function(nodes, edges, dimNames, onlyNovel = F){
  #fix values of confidence levels
  nodes$Confidence = conf.f(nodes$Confidence)
  
  # change names of from and to to name and imports
  colnames(edges)[1:2] = c('name', 'imports')
  
  # creates a larger list designating node connection to all other nodes
  group1 = lapply(split(edges, edges$name), `[`, 2:ncol(edges))
  group2 = lapply(split(edges, edges$imports), `[`, c(1,3:ncol(edges)))
  group2 = lapply(group2, setNames, c('imports', 'weights', 'datasource'))
  L = merge_lists(group1, group2)
  
  # orders pathways for 1st and 2nd degreee
  if(!onlyNovel){
    L = orderL(L, dimNames)
  }
  
  # json list formation where each entry of the list is a dataframe of connections for a node
  df_L = list()
  for(i in 1:length(L)){df_L = append(df_L, list(change_form(L, i, nodes)))}
  return(df_L)
}





