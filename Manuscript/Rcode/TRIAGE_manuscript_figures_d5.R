##################################
## Figures for TRIAGE manuscript
##################################

library(hexbin)
library(tidyr)
library(dplyr)
library(ggplot2)
library(xlsx)
library(VennDiagram)
library(org.Hs.eg.db)
library(ggthemes)
library(hexbin)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(scales)
library(circlize)



## Basics
#Figures working directory
figures.wd <- "~/OneDrive\ -\ National\ Institutes\ of\ Health/TRIAGE/Manuscript/Text/Drafts/Full_Draft_d4/R_figures/"
results.wd <- "~/OneDrive\ -\ National\ Institutes\ of\ Health/TRIAGE/Manuscript/Text/Drafts/Full_Draft_d4/R_figures/"
######## HIV Reported Hits -----------
#Set colors

postval.color <- "#A9C353"
highscore.color <- "#8A364C"


HIV.wd <- "/Users/sakatz/Documents/Analysis/HIV/TRIAGEHIV_May2019/TRIAGE/TRIAGEinput/Files/"
setwd(HIV.wd)
#Get Three HIV TRIAGE input files
screen.names <- c("zhou", "brass", "konig")
for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  assign(paste0(screen.name, ".df"), read.csv(paste0(screen.name, ".triageinput.csv"), stringsAsFactors = F))
  assign(paste0(screen.name, ".validated.matrix"), as.matrix(get(paste0(screen.name, ".df"))[which(get(paste0(screen.name, ".df"))[["validated.triage"]] == 1), c("EntrezID")]))
}
#Get ven of reported hits
fill.color <- c(postval.color, postval.color, postval.color)
setwd(paste0(figures.wd))
tophits_graph <- venn.diagram(
  x = list(get(paste0("zhou.validated.matrix")), get(paste0("brass.validated.matrix")), get(paste0("konig.validated.matrix"))),
  category.names = c("Zhou" , "Brass" , "Konig"),
  filename = paste0("HIV_postval", "_venn.jpeg"),
  output = TRUE ,
  imagetype="png" ,
  height = 480 ,
  width = 480 ,
  resolution = 300,
  compression = "lzw",
  lwd = 0.2,
  lty = 'solid',
  col = 'black',
  fill = fill.color,
  alpha = 0.5,
  label.col = c("#626262", "black", "#626262", "black", "black", "black", "#626262"),
  cex = 1,
  fontface = c("plain", "bold", "plain", "bold","bold","bold", "plain"),
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


#libraries
triage.gen.wd <- "~/Documents/Analysis/HIV/TRIAGEHIV_May2019"
hiv.raw.wd <- paste0(triage.gen.wd, "/OriginalData/AuthorProvided") 

# get files ---> for each screen get name.raw.df, name.bhaskar.df, name.card.df

#Zhou
zhou.raw.df <- read.csv(paste0(hiv.raw.wd, "/Zhou/MerckDataOriginal_genes.csv"), stringsAsFactors = F)

#Brass
brass.raw.df <- read.csv(paste0(hiv.raw.wd, "/Brass/SMARTpoolHIVscreenforIainFeb26.csv"), stringsAsFactors = F)

#Konig
konig.raw.df <- read.csv(paste0(hiv.raw.wd, "/Konig/ChandaOriginalData.csv"), stringsAsFactors = F)

#####- Normalization
#Zhou-
#Log transfre inhibition to expression levels and Zscore for viability and remove non human genes
#Define humna gene list
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)

zhou.normalized.df <-  zhou.raw.df %>%
  mutate(X48Hexpression = log(1 - (0.01 * X48HpercINH))) %>%
  mutate(X96Hexpression = log(1 - (0.01 * X96HpercINH))) %>%
  subset(EntrezID %in% mapped_genes) %>%
  mutate(X48hViab.Zscore = as.numeric(scale(percViability48H, center = TRUE, scale = TRUE))) %>%
  mutate(X96hViab.Zscore = as.numeric(scale(percViability96H, center = TRUE, scale = TRUE)))

#colnames
zhou.raw.name <- "X48HpercINH"
zhou.normalized.name <- "X48Hexpression"
zhou.confidence.name <- "X48hViab.Zscore"
zhou.x.ticks <- c(-3, -2,-1, 0, 1, 2, 3)
zhou.y.ticks <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)

#Brass-

# - Remove rows with no Entrez, log normalize cell count.
brass.normalized.df <- brass.raw.df %>%
  subset(is.na(EntrezID) == F) %>%
  subset(EntrezID != 0) %>%
  mutate(Log10CellNumber = log10(NormalizedCellNumber)) %>%
  mutate(UniqueID = paste0(Plate, Well)) #Unique identifier



#-nromalize percent infected and cell count by plate 
plates <- unique(brass.normalized.df$Plate)   #Get list of plate ID

brass.normalized.df$PercInfected.Zscore <- NA #Name column for normalized values

PlateNormCol <- which(colnames(brass.normalized.df) == "PercInfected.Zscore")
DataCol <- which(colnames(brass.normalized.df) == "NormalizedPercentInfected")

for(i in 1:length(plates)){
  indPlate <- which(brass.normalized.df$Plate == plates[i])
  
  brass.normalized.df[indPlate, PlateNormCol] <- (brass.normalized.df[indPlate,DataCol] -  mean(brass.normalized.df[indPlate,DataCol],na.rm = TRUE))/
    sd(brass.normalized.df[indPlate,DataCol],na.rm = TRUE)
}

brass.normalized.df$CellNumber.Zscore <- NA

PlateNormCol <- which(colnames(brass.normalized.df) == "CellNumber.Zscore")
DataCol <- which(colnames(brass.normalized.df) == "Log10CellNumber")

for(i in 1:length(plates)){
  indPlate <- which(brass.normalized.df$Plate == plates[i])
  
  brass.normalized.df[indPlate, PlateNormCol] <- (brass.normalized.df[indPlate,DataCol] -  mean(brass.normalized.df[indPlate,DataCol],na.rm = TRUE))/
    sd(brass.normalized.df[indPlate,DataCol],na.rm = TRUE)
}


brass.raw.name <- "NormalizedPercentInfected"
brass.normalized.name <- "PercInfected.Zscore"
brass.confidence.name <- "CellNumber.Zscore"
brass.x.ticks <- c(-5, -4, -3, -2,-1, 0, 1, 2, 3, 4)
brass.y.ticks <- c(-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3)

#- Konig
konig.normalized.df <- konig.raw.df %>%
  tidyr::separate(Well_ID, c("Plate", "Well") ) %>%
  mutate(PlateID = paste(Plate, substring(Library, 4), sep = "")) %>%
  mutate(negLogP = LogP*-1) %>%
  mutate(LogScore = log(Score)) %>%
  mutate(UniqueID = paste0(PlateID, "_", Well)) #Unique identifier


konig.normalized.df$Plate <- NULL  

names(konig.normalized.df)[which(colnames(konig.normalized.df) == "Gene_ID")] <- "EntrezID"

#-nromalize percent infected and cell count by plate 
plates <- unique(konig.normalized.df$PlateID)   #Get list of plate ID

konig.normalized.df$Zscore <- NA #Name column for normalized values

PlateNormCol <- which(colnames(konig.normalized.df) == "Zscore")
DataCol <- which(colnames(konig.normalized.df) == "LogScore")

for(i in 1:length(plates)){
  indPlate <- which(konig.normalized.df$PlateID == plates[i])
  
  konig.normalized.df[indPlate, PlateNormCol] <- (konig.normalized.df[indPlate,DataCol] -  mean(konig.normalized.df[indPlate,DataCol],na.rm = TRUE))/
    sd(konig.normalized.df[indPlate,DataCol],na.rm = TRUE)
}

konig.raw.name <- "Score"
konig.normalized.name <- "Zscore"
konig.confidence.name <- "negLogP"
konig.x.ticks <- c(-5, -4, -3, -2,-1, 0, 1, 2, 3, 4, 5)
konig.y.ticks <- c( -1, 0, 1, 2, 3, 4, 5, 6)


####### Figures
screen.names <- c("zhou", "brass", "konig")
raw.columns <- c(zhou.raw.name, brass.raw.name, konig.raw.name)
normalized.columns <- c(zhou.normalized.name, brass.normalized.name, konig.normalized.name)
confidence.columns <- c(zhou.confidence.name, brass.confidence.name, konig.confidence.name)
confidence.names <- c("Viability", "Viability", "-LogP")




# count versus score
for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  
  assign(paste0(screen.name, ".plot"), ggplot(get(paste0(screen.name, ".normalized.df")))+
           geom_histogram(aes(x = get(normalized.columns[s])), binwidth = .1, 
                          col="paleturquoise1", 
                          size=.1, fill="black")+
           ylab("Frequency") + xlab("Zscore")+
           ggtitle(paste0(toupper(screen.name), ": Score Distribution"))+
           theme(legend.position="none", axis.text.y = element_text(size=15),
                 axis.text.x = element_text(size=10),
                 axis.title=element_text(size=20,face="bold"),
                 axis.line.y = element_line(colour = 'black', size = 2),
                 axis.line.x = element_line(colour = 'black', size = 2),
                 panel.background = element_blank()))
  
  setwd(paste0(figures.wd))
  ggsave(paste0(screen.name, ".scorefrequency.png"), get(paste0(screen.name, ".plot")))
  
}


#- get median zscore and median confidence measure for each screen
for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  #aggregate the score
  agg.score <- aggregate(get(normalized.columns[s]) ~ EntrezID, get(paste0(screen.name, ".normalized.df")), median)
  names(agg.score)[2] <- paste0((normalized.columns[s]))
  agg.conf <- aggregate(get(confidence.columns[s]) ~ EntrezID, get(paste0(screen.name, ".normalized.df")), median)
  names(agg.conf)[2] <- paste0((confidence.columns[s]))
  
  assign(paste0(screen.name, ".median.df"), merge(agg.score, agg.conf))
  
  
}





#Get list of hits for each study
Supplement.WD <- "/Users/sakatz/Documents/Analysis/HIV/SupplemntaryTables/"

###  Get Author Selected Hits and create matrices 
# Number of hits per screen
zhou.hits.count <- 232
brass.hits.count <- 281
konig.hits.count <- 295

setwd(Supplement.WD)
AuthorSelectedHits.df <- read.csv("AuthorSelectedHits_HIV.csv", stringsAsFactors = F)

zhou.top.hits <- matrix(AuthorSelectedHits.df$ZhouHits[1:zhou.hits.count])
brass.top.hits <- matrix(AuthorSelectedHits.df$BrassHits[1:brass.hits.count])
konig.top.hits <- matrix(AuthorSelectedHits.df$KonigHits[1:konig.hits.count])



screen.confidence.thresholds <- c(-2.00, -2.00, 1.3)
#- Defines sizes of subsets
HighConf.size <- 400
MedConf.size <- 1000

#- get cutoff for 
for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  
  #- subset dataframe 
  #-- Ranking based on CARD and confidence measure
  subset.df <- filter(get(paste0(screen.name, ".median.df")), get(confidence.columns[s]) > screen.confidence.thresholds[s]) %>%
    arrange(get(normalized.columns[s]))
  
  assign(paste0(screen.name, ".scorecutoff"), subset.df[[normalized.columns[s]]][HighConf.size])
}

# Visualize



#with top hits highlight
# Paremeters
shape = 16
size = 5
for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  
  assign(paste0(screen.name, ".plot"), ggplot(get(paste0(screen.name, ".median.df")))+
           geom_rect(aes(xmin=-Inf,xmax=get(paste0(screen.name, ".scorecutoff")),ymin=screen.confidence.thresholds[s],ymax=Inf), alpha=1,fill = highscore.color)+
           geom_point(aes(x = get(normalized.columns[s]), y = get(confidence.columns[s]), colour = EntrezID %in% get(paste0(screen.name, ".top.hits"))), shape = shape, size = size)+
           geom_point(data = subset(get(paste0(screen.name, ".median.df")), EntrezID %in% get(paste0(screen.name, ".top.hits"))),
                      aes(x = get(normalized.columns[s]), y =get(confidence.columns[s]), colour = EntrezID %in% get(paste0(screen.name, ".top.hits"))), shape = shape, size = size)+
           scale_colour_manual(name = 'Selected Hits', values = setNames(c(postval.color,'black'),c(T, F))) +
           scale_shape_manual()+
           scale_x_continuous(breaks = get(paste0(screen.name, ".x.ticks")))+
           scale_y_continuous(breaks = get(paste0(screen.name, ".y.ticks")))+
           xlab('Zscore') + ylab(confidence.names[s])+
           ggtitle(paste0(toupper(screen.name), ": Median Score & Confidence"))+
           theme(legend.position="none", axis.text.y = element_text(size=15),
                 axis.text.x = element_text(size=10),
                 axis.title=element_text(size=20,face="bold"),
                 axis.line.y = element_line(colour = 'black', size = 2),
                 axis.line.x = element_line(colour = 'black', size = 2),
                 panel.background = element_blank(),
                 aspect.ratio = 1))
  
  setwd(paste0(figures.wd))
  
  ggsave(paste0(screen.name, ".postval.jpeg"), get(paste0(screen.name, ".plot")))
}


# Add hits density factors.
# Paremeters
shape = 16
size = 5
for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  
  assign(paste0(screen.name, ".plot"), ggplot(get(paste0(screen.name, ".median.df")))+
           geom_rect(aes(xmin=-Inf,xmax=get(paste0(screen.name, ".scorecutoff")),ymin=screen.confidence.thresholds[s],ymax=Inf), alpha=1,fill = "red")+
           geom_hex(aes(x = get(normalized.columns[s]), y = get(confidence.columns[s]), colour = EntrezID %in% get(paste0(screen.name, ".top.hits"))), bins = 70)+
           # geom_hex(data = subset(get(paste0(screen.name, ".median.df")), EntrezID %in% get(paste0(screen.name, ".top.hits"))),
           #            aes(x = get(normalized.columns[s]), y =get(confidence.columns[s]), colour = EntrezID %in% get(paste0(screen.name, ".top.hits"))), bins = 70)+
           scale_fill_continuous(type = "gradient", name = "Density") +
           #scale_fill_distiller() +
           scale_colour_manual(name = 'Selected Hits', values = setNames(c('green','black'),c(T, F)), guide = F) +
           scale_shape_manual()+
           xlab('Zscore') + ylab(confidence.names[s])+
           ggtitle(paste0(toupper(screen.name), ": Median Score & Confidence"))+
           theme(legend.justification ="top", legend.key.height = unit(2, "cm"),
                 legend.key.width = unit(0.7, "cm"),
                 legend.title = element_text(size=30, face="bold"),
                 legend.text = element_text(size=25),
                 axis.text.y = element_text(size=15),
                 axis.text.x = element_text(size=10),
                 axis.title=element_text(size=20,face="bold"),
                 axis.line.y = element_line(colour = 'black', size = 2),
                 axis.line.x = element_line(colour = 'black', size = 2),
                 panel.background = element_blank()))
  
  setwd(paste0(figures.wd))
  
  ggsave(paste0(screen.name, ".postval.density.jpeg"), get(paste0(screen.name, ".plot")))
}



# Venn diagram of overlap
#Get Three HIV TRIAGE hit matrices files

screen.names <- c("zhou", "brass", "konig")

for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  
  assign(paste0(screen.name, ".normalized.matrix"), as.matrix(get(paste0(screen.name, ".df"))[which(get(paste0(screen.name, ".df"))[["normalized.triage"]] == 1), c("EntrezID")]))
}


#Get ven of reported hits

fill.color <- c(highscore.color, highscore.color, highscore.color)



setwd(paste0(figures.wd))

tophits_graph <- venn.diagram(
  x = list(get(paste0("zhou.normalized.matrix")), get(paste0("brass.normalized.matrix")), get(paste0("konig.normalized.matrix"))),
  category.names = c("Zhou" , "Brass" , "Konig"),
  filename = paste0("HIV_topscore", "_venn.jpeg"),
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 0.2,
  lty = 'solid',
  col = 'black',
  fill = fill.color,
  alpha = 0.5,
  label.col = c("#626262", "black", "#626262", "black", "black", "black", "#626262"),
  cex = 1,
  fontface = c("plain", "bold", "plain", "bold","bold","bold", "plain"),
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)



###Make Bar Chart of Pval and Shared hits

#get overlap results
source("~/Documents/Functions/TwoScreenOverlap.R")

Zhou_Brass.val.overlap <- OVERLAP.function(zhou.df, brass.df, "EntrezID", "validated.triage", 1, "equal")
Brass_Konig.val.overlap <- OVERLAP.function(brass.df, konig.df, "EntrezID", "validated.triage", 1, "equal")
Zhou_Konig.val.overlap <- OVERLAP.function(zhou.df, konig.df, "EntrezID", "validated.triage", 1, "equal")

Zhou_Brass.top.overlap <- OVERLAP.function(zhou.df, brass.df, "EntrezID", "normalized.triage", 1, "equal")
Brass_Konig.top.overlap <- OVERLAP.function(brass.df, konig.df, "EntrezID", "normalized.triage", 1, "equal")
Zhou_Konig.top.overlap <- OVERLAP.function(zhou.df, konig.df, "EntrezID", "normalized.triage", 1, "equal")


### Get hit selection set size
Zhou_top_set <- length(zhou.df$EntrezID[which(zhou.df$normalized.triage == 1)])
Brass_top_set <- length(brass.df$EntrezID[which(brass.df$normalized.triage == 1)])
Konig_top_set <- length(konig.df$EntrezID[which(konig.df$normalized.triage == 1)])

Zhou_val_set <- length(zhou.df$EntrezID[which(zhou.df$validated.triage == 1)])
Brass_val_set <- length(brass.df$EntrezID[which(brass.df$validated.triage == 1)])
Konig_val_set <- length(konig.df$EntrezID[which(konig.df$validated.triage == 1)])

#### Figures
# write dataframe
pval.df <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                            , "HitSelection" = "High Scoring Hits"
                            , "LogpVal" = c(log10(Zhou_Brass.top.overlap[[3]])*-1
                                            ,log10(Brass_Konig.top.overlap[[3]])*-1
                                            , log10(Zhou_Konig.top.overlap[[3]])*-1
                            ))
                 ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                             , "HitSelection" = "Post-Validation Hits"
                             , "LogpVal" = c(log10(Zhou_Brass.val.overlap[[3]])*-1
                                             ,log10(Brass_Konig.val.overlap[[3]])*-1
                                             , log10(Zhou_Konig.val.overlap[[3]])*-1
                             )))


# Faceting
pval.figure <- ggplot(pval.df, aes(y=LogpVal, x=Comparison, color=HitSelection, fill=HitSelection)) + 
  geom_bar( stat="identity", width = 0.7) +  
  #scale_color_manual(values = c("#FF7676", "#05FDFF"))+
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = c(highscore.color, postval.color))+
  facet_wrap(~HitSelection) +
  scale_y_continuous(limits=c(0, 5), breaks = c(0, log10(0.05)*-1, 3, 4, 5),
                     labels = c(round(0, 1), round(log10(0.05)*-1, 2), c(round(c(3,4,5), 1)))) +
  geom_hline(yintercept=log10(0.05)*-1, linetype="dashed", 
             color = "black", size=0.5)+
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=12, colour = "black", angle = 45, hjust = 1),
        axis.title=element_text(size=15,face="bold"),
        axis.text.y.left = element_text(size=15),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        aspect.ratio = 2.5/1,
        panel.background = element_blank())

setwd(paste0(figures.wd))

ggsave(paste0("pval_comparison_hiv.png"), pval.figure)


#Share hits graph
# write dataframe
sharedhits.df <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                  , "HitSelection" = "High Scoring Hits"
                                  , "SharedHits" = c(Zhou_Brass.top.overlap[[2]]
                                                     ,Brass_Konig.top.overlap[[2]]
                                                     ,Zhou_Konig.top.overlap[[2]]
                                  ))
                       ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                   , "HitSelection" = "Post-Validation Hits"
                                   , "SharedHits" = c(Zhou_Brass.val.overlap[[2]]
                                                      ,Brass_Konig.val.overlap[[2]]
                                                      ,Zhou_Konig.val.overlap[[2]]
                                   )))


# Faceting
sharedhits.figure <- ggplot(sharedhits.df, aes(y=SharedHits, x=Comparison, color=HitSelection, fill=HitSelection)) + 
  geom_bar( stat="identity", width = 0.7) +   
  #scale_color_manual(values = c("#FF7676", "#05FDFF"))+
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = c(highscore.color, postval.color))+
  facet_wrap(~HitSelection) +
  # scale_y_continuous(limits=c(0, 10), breaks = c(0, log10(0.05)*-1, 2, 4, 6, 8, 10),
  #                     labels = c(round(0, 1), round(log10(0.05)*-1, 2), c(round(c(2, 4, 6, 8, 10), 1)))) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=12, colour = "black", angle = 45, hjust = 1),
        axis.title=element_text(size=15,face="bold"),
        axis.text.y.left = element_text(size=15),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        aspect.ratio = 2.5/1,
        panel.background = element_blank())

setwd(paste0(figures.wd))

ggsave(paste0("sharedhits_hiv.png"), sharedhits.figure)


#Set size graph
# write dataframe
setsize.df <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                               , "HitSelection" = "High Scoring Hits"
                               , "AverageSetSize" = c(mean(Zhou_top_set, Brass_top_set)
                                                      ,mean(Brass_top_set, Konig_top_set)
                                                      ,mean(Konig_top_set, Zhou_top_set)
                               ))
                    ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                , "HitSelection" = "Post-Validation Hits"
                                , "AverageSetSize" = c(mean(Zhou_val_set, Brass_val_set)
                                                       ,mean(Brass_val_set, Konig_val_set)
                                                       ,mean(Konig_val_set, Zhou_val_set)
                                )))


# Faceting
setsize.figure <- ggplot(setsize.df, aes(y=AverageSetSize, x=Comparison, color=HitSelection, fill=HitSelection)) + 
  geom_bar( stat="identity", width = 0.7) +  
  #scale_color_manual(values = c("#FF7676", "#05FDFF"))+
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = c(highscore.color, postval.color))+
  facet_wrap(~HitSelection) +
  # scale_y_continuous(limits=c(0, 10), breaks = c(0, log10(0.05)*-1, 2, 4, 6, 8, 10),
  #                     labels = c(round(0, 1), round(log10(0.05)*-1, 2), c(round(c(2, 4, 6, 8, 10), 1)))) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=12, colour = "black", angle = 45, hjust = 1),
        axis.title=element_text(size=15,face="bold"),
        axis.text.y.left = element_text(size=15),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        aspect.ratio = 2.5/1,
        panel.background = element_blank())

setwd(paste0(figures.wd))

ggsave(paste0("setsize_hiv.png"), setsize.figure)

############## FIGURE 2 #############
for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  
  assign(paste0(screen.name, ".plot"), ggplot(get(paste0(screen.name, ".median.df")))+
           geom_rect(aes(xmin=-Inf,xmax=get(paste0(screen.name, ".scorecutoff")),ymin=screen.confidence.thresholds[s],ymax=Inf), alpha=1,fill = highscore.color)+
           geom_point(aes(x = get(normalized.columns[s]), y = get(confidence.columns[s]), colour = EntrezID %in% get(paste0(screen.name, ".top.hits"))), shape = shape, size = size)+
           geom_point(data = subset(get(paste0(screen.name, ".median.df")), EntrezID %in% get(paste0(screen.name, ".top.hits"))),
                      aes(x = get(normalized.columns[s]), y =get(confidence.columns[s]), colour = EntrezID %in% get(paste0(screen.name, ".top.hits"))), shape = shape, size = size)+
           scale_colour_manual(name = 'Selected Hits', values = setNames(c(postval.color,'black'),c(T, F))) +
           scale_shape_manual()+
           scale_x_continuous(breaks = get(paste0(screen.name, ".x.ticks")))+
           scale_y_continuous(breaks = get(paste0(screen.name, ".y.ticks")))+
           xlab('Zscore') + ylab(confidence.names[s])+
           ggtitle(paste0(toupper(screen.name), ": Median Score & Confidence"))+
           theme(legend.position="none", axis.text.y = element_text(size=15),
                 axis.text.x = element_text(size=10),
                 axis.title=element_text(size=20,face="bold"),
                 axis.line.y = element_line(colour = 'black', size = 2),
                 axis.line.x = element_line(colour = 'black', size = 2),
                 panel.background = element_blank(),
                 aspect.ratio = 1))
  
  setwd(paste0(figures.wd))
  
  ggsave(paste0(screen.name, ".postval.jpeg"), get(paste0(screen.name, ".plot")))
}



tier <- "normalized"

#Set colors
HC.color <- "#DE8731"
MC.color <- "#7344b7" 
LC.color <- "black"

tiers <- c("normalized", "validated")

for (t in 1:length(tiers)){
  tier <- tiers[t]
  
  for (s in 1:length(screen.names)){
    screen.name <- screen.names[s]
    
    HiConf.hits <- as.matrix(get(paste0(screen.name, ".df"))[which(get(paste0(screen.name, ".df"))[[paste0(tier, ".triage")]] == 1), c("EntrezID")])
    MedConf.hits <- as.matrix(get(paste0(screen.name, ".df"))[which(get(paste0(screen.name, ".df"))[[paste0(tier, ".triage")]] == 0.5), c("EntrezID")])
    
    
    assign(paste0(screen.name, ".plot"), ggplot(get(paste0(screen.name, ".df")))+
             geom_point(aes(x = get(normalized.columns[s]), y = get(confidence.columns[s]), colour = as.character(get(paste0(tier, ".triage")))), shape = shape, size=size)+
             geom_point(data = subset(get(paste0(screen.name, ".df")), get(paste0(tier, ".triage")) == 0.5),
                        aes(x = get(normalized.columns[s]), y = get(confidence.columns[s]), colour = as.character(get(paste0(tier, ".triage")))), shape = shape, size=size)+
             geom_point(data = subset(get(paste0(screen.name, ".df")), get(paste0(tier, ".triage")) == 1),
                        aes(x = get(normalized.columns[s]), y = get(confidence.columns[s]), colour = as.character(get(paste0(tier, ".triage")))), shape = shape, size=size)+
             #scale_colour_manual(name = 'Selected Hits', values = setNames(c('red', 'blue','black')), c("High Conf", "MidConf", "LowConf")) +
             scale_colour_manual(name = 'Selected Hits', values = c("1" = HC.color, "0.5" = MC.color, "0" = LC.color), guide = F) +
             scale_shape_manual()+
             scale_x_continuous(breaks = get(paste0(screen.name, ".x.ticks")))+
             scale_y_continuous(breaks = get(paste0(screen.name, ".y.ticks")))+
             xlab('Zscore') + ylab(confidence.names[s])+
             ggtitle(paste0(toupper(screen.name), ": post-", tier, " TRIAGE input"))+
             theme(legend.position="none", axis.text.y = element_text(size=15),
                   axis.text.x = element_text(size=10),
                   axis.title=element_text(size=20,face="bold"),
                   axis.line.y = element_line(colour = 'black', size = 2),
                   axis.line.x = element_line(colour = 'black', size = 2),
                   panel.background = element_blank(),
                   aspect.ratio = 1))
    
    
    setwd(paste0(results.wd, "SegmentationFigures/"))
    ggsave(paste0(screen.name, ".", tier, "triageinput.jpeg"), get(paste0(screen.name, ".plot")))
  }
}
  







### Get csv with with overlap and permutation test data

comparison.analysis.df <- read.csv(paste0("~/Documents/Analysis/HIV/TRIAGEHIV_May2019/Comparison/PermutationTests/ComparisonAnalysis_1000permTest.csv"), stringsAsFactors = FALSE)

###-- Seperate Dataframes by analysis and filter by analysis you want to include in figure

analysis.included <- c("UGP", "HC", "HC_PA", "HC_PA2t", "HC_NA", "SERIAL_PA2t", "TRIAGE_PA2t")

overlap.df <- comparison.analysis.df[which(grepl("Overlap", comparison.analysis.df$X) == T), ]
overlap.df <- comparison.analysis.df[which(comparison.analysis.df$X %in% paste0(analysis.included, "_Overlap")), ]

pval.df <- comparison.analysis.df[which(grepl("negLogP", comparison.analysis.df$X) == T), ]
pval.df <- comparison.analysis.df[which(comparison.analysis.df$X %in% paste0(analysis.included, "_negLogP")), ]

perm_test.df <- comparison.analysis.df[which(grepl("Random.perc", comparison.analysis.df$X) == T), ]
perm_test.df <- comparison.analysis.df[which(comparison.analysis.df$X %in% paste0(analysis.included, "_Random.perc")), ]

#- Filter for the analysis you want to include 
##- Find max for each dataframe
overlap.max <- (max(overlap.df[, -1]) + 1)
pval.max <- ceiling(max(pval.df[, c(-1, -5)]))

overlap.yticks <- c(ceiling((overlap.max/3)*1), ceiling((overlap.max/3)*2), overlap.max)
pval.yticks <- c(ceiling((pval.max/3)*1), ceiling((pval.max/3)*2), pval.max)



###### Post Validaion Hits -----

#Create data frame for graph - pVal and shared hits
analysis.type <- "UGP"

postval.hits.pval <- data.frame("Comparison" = c("Zhou + Brass", "Brass + Konig", "Zhou + Konig")
                                , "HitSelection" = "High Scoring Hits"
                                , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                ))

postval.hits.shared <- data.frame("Comparison" = c("Zhou + Brass", "Brass + Konig", "Zhou + Konig")
                                  , "HitSelection" = "High Scoring Hits"
                                  , "SharedHits" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                  ))



## Create figure

levels(postval.hits.pval$Comparison) <- gsub(" ", "\n", levels(postval.hits.pval$Comparison))

pval.figure <- ggplot(postval.hits.pval, aes(y=LogpVal, x=Comparison)) + 
  geom_bar( stat="identity", width = 0.7, color= "darkred", fill="darkred") +    
  #facet_wrap(~HitSelection) +
  ylab("-Log p Value")+
  scale_y_continuous(limits=c(0, pval.max), breaks = c(0, pval.yticks),
                     labels = c(0, pval.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=35),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        aspect.ratio = 1.5/1,
        panel.background = element_blank())

print(pval.figure)

setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("test_pval_comparison_postval.png"), pval.figure)

levels(postval.hits.shared$Comparison) <- gsub(" ", "\n", levels(postval.hits.shared$Comparison))

overlap.figure <- ggplot(postval.hits.shared, aes(y=SharedHits, x=Comparison)) + 
  geom_bar( stat="identity", width = 0.7, color= "darkblue", fill="darkblue") +    
  #facet_wrap(~HitSelection) +
  ylab("Shared Hits")+
  scale_y_continuous(limits=c(0, overlap.max), breaks = c(0, overlap.yticks),
                     labels = c(0, overlap.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=35),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        aspect.ratio = 1.5/1,
        panel.background = element_blank())

setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("overlap_comparison_postval.png"), overlap.figure)


###### High Scoring Hits -----

#Create data frame for graph - pVal and shared hits
analysis.type <- "HC"

highscore.hits.pval <- data.frame("Comparison" = c("Zhou + Brass", "Brass + Konig", "Zhou + Konig")
                                  , "HitSelection" = "High Scoring Hits"
                                  , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                  ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                  ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                  ))

highscore.hits.shared <- data.frame("Comparison" = c("Zhou + Brass", "Brass + Konig", "Zhou + Konig")
                                    , "HitSelection" = "High Scoring Hits"
                                    , "SharedHits" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                       ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                       ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                    ))



## Create figure
levels(highscore.hits.pval$Comparison) <- gsub(" ", "\n", levels(highscore.hits.pval$Comparison))

pval.figure <- ggplot(highscore.hits.pval, aes(y=LogpVal, x=Comparison)) + 
  geom_bar( stat="identity", width = 0.7, color= "darkred", fill="darkred") +    
  #facet_wrap(~HitSelection) +
  ylab("-Log p Value")+
  scale_y_continuous(limits=c(0, pval.max), breaks = c(0, pval.yticks),
                     labels = c(0, pval.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=35),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        aspect.ratio = 1.5/1,
        panel.background = element_blank())

setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("pval_comparison_highscore.png"), pval.figure)

levels(highscore.hits.shared$Comparison) <- gsub(" ", "\n", levels(highscore.hits.shared$Comparison))

overlap.figure <- ggplot(highscore.hits.shared, aes(y=SharedHits, x=Comparison)) + 
  geom_bar( stat="identity", width = 0.7, color= "darkblue", fill="darkblue") +    
  #facet_wrap(~HitSelection) +
  ylab("Shared Hits")+
  scale_y_continuous(limits=c(0, overlap.max), breaks = c(0, overlap.yticks),
                     labels = c(0, overlap.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=35),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        aspect.ratio = 1.5/1,
        panel.background = element_blank())

setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("overlap_comparison_highscore.png"), overlap.figure)

###### Set plot colors --------
### Set color
overlap.color <- "#9ED867"
pval.color <- "#3B82A2"


pval.plot.colors <- c("grey", "black", pval.color)
overlap.plot.colors <- c("grey", "black", overlap.color)

pval.plot.colors2 <- c("grey", "black", pval.color, "#4d7689")
overlap.plot.colors2 <- c("grey", "black", overlap.color, "#9fbf82")

pval.border <- c("black", "black", "black")
overlap.border <- c("black", "black", "black")

pval.border2 <- c("black", "black", "black", "black")
overlap.border2 <- c("black", "black", "black", "black")

border.size <- 1

###### Two tier pathway analysis Hits -----

#Create data frame for graph - pVal and shared hits
analysis.type <- "HC_PA"
analysis.name <- "Pathway Analysis"
analysis.file.name <- "pathway"

#Create levels for factors
Levels <- c("High Scoring Hits", "Post Validation Hits", analysis.name)

hits.pval <- rbind(data.frame("Comparison" = c("Zhou + Brass", "Brass + Konig", "Zhou + Konig")
                              , "HitSelection" = analysis.name
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              )),
                   data.frame("Comparison" = c("Zhou + Brass", "Brass + Konig", "Zhou + Konig")
                              , "HitSelection" = "High Scoring Hits"
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              ))
                   ,data.frame("Comparison" = c("Zhou + Brass", "Brass + Konig", "Zhou + Konig")
                               , "HitSelection" = "Post Validation Hits"
                               , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                               )))
#Introduce linebreaks
levels(hits.pval$Comparison) <- gsub(" ", "\n", levels(hits.pval$Comparison))
#Order factor of groups
hits.pval$HitSelection <- factor(hits.pval$HitSelection, levels = Levels)

#Shared hits wih pval

hits.shared <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                , "HitSelection" = analysis.name
                                , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "High Scoring Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "Post Validation Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 )))
#Introduce linebreaks
levels(hits.shared$Comparison) <- gsub(" ", "\n", levels(hits.shared$Comparison))
levels(hits.shared$Comparison) <- gsub("-", "+", levels(hits.shared$Comparison))

hits.shared$HitSelection <- factor(hits.shared$HitSelection, levels = Levels)

#Assign star designation for significance
hits.shared$star <- ""
hits.shared$star[hits.shared$Significance > .05 & hits.shared$HitSelection == analysis.name]  <- "ns"
hits.shared$star[hits.shared$Significance <= .05 & hits.shared$HitSelection == analysis.name]  <- "*"
hits.shared$star[hits.shared$Significance <= .01 & hits.shared$HitSelection == analysis.name]  <- "**"




## Create figure



pval.figure <- ggplot(hits.pval, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection) +
  ylab("-Log p Value")+
  scale_color_manual(values = pval.border)+
  scale_fill_manual(breaks = Levels,
                    values = pval.plot.colors)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  scale_y_continuous(limits=c(0, pval.max), breaks = c(0, pval.yticks),
                     labels = c(0, pval.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())
print(pval.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))


ggsave(paste0("pval_comparison_", analysis.file.name, ".png"), pval.figure)


overlap.figure <- ggplot(hits.shared, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection) +
  ylab("Shared Hits")+
  scale_color_manual(values = overlap.border)+
  scale_fill_manual(breaks = Levels,
                    values = overlap.plot.colors)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  geom_text(aes(label=star), colour = "black", position = position_dodge(width = 0.85), vjust=0, size=10) +
  scale_y_continuous(limits=c(0, overlap.max), breaks = c(0, overlap.yticks),
                     labels = c(0, overlap.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())

print(overlap.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))


ggsave(paste0("overlap_comparison_", analysis.file.name, ".png"), overlap.figure)

###### Three tier pathway analysis Hits -----

#Create data frame for graph - pVal and shared hits
analysis.type <- "HC_PA2t"
analysis.name <- "Pathway Analysis"
analysis.file.name <- "pathway2t"


#Create levels for factors
Levels <- c("High Scoring Hits", "Post Validation Hits", analysis.name)

hits.pval <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = analysis.name
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              )),
                   data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = "High Scoring Hits"
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              ))
                   ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                               , "HitSelection" = "Post Validation Hits"
                               , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                               )))
#Introduce linebreaks
levels(hits.pval$Comparison) <- gsub(" ", "\n", levels(hits.pval$Comparison))
levels(hits.pval$Comparison) <- gsub("-", "+", levels(hits.pval$Comparison))
#Order factor of groups
hits.pval$HitSelection <- factor(hits.pval$HitSelection, levels = Levels)

#Shared hits wih pval

hits.shared <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                , "HitSelection" = analysis.name
                                , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "High Scoring Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "Post Validation Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 )))

#Introduce linebreaks
levels(hits.shared$Comparison) <- gsub(" ", "\n", levels(hits.shared$Comparison))
levels(hits.shared$Comparison) <- gsub("-", "+", levels(hits.shared$Comparison))

hits.shared$HitSelection <- factor(hits.shared$HitSelection, levels = Levels)

#Assign star designation for significance
hits.shared$star <- ""
hits.shared$star[hits.shared$Significance > .05 & hits.shared$HitSelection == analysis.name]  <- "ns"
hits.shared$star[hits.shared$Significance <= .05 & hits.shared$HitSelection == analysis.name]  <- "*"
hits.shared$star[hits.shared$Significance <= .01 & hits.shared$HitSelection == analysis.name]  <- "**"




## Create figure



pval.figure <- ggplot(hits.pval, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection) +
  ylab("-Log p Value")+
  scale_color_manual(values = pval.border)+
  scale_fill_manual(breaks = Levels,
                    values = pval.plot.colors)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  scale_y_continuous(limits=c(0, pval.max), breaks = c(0, pval.yticks),
                     labels = c(0, pval.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())
print(pval.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("pval_comparison_", analysis.file.name, ".png"), pval.figure)


overlap.figure <- ggplot(hits.shared, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection) +
  ylab("Shared Hits")+
  scale_color_manual(values = overlap.border)+
  scale_fill_manual(breaks = Levels,
                    values = overlap.plot.colors)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  geom_text(aes(label=star), colour = "black", position = position_dodge(width = 0.85), vjust=0, size=10) +
  scale_y_continuous(limits=c(0, overlap.max), breaks = c(0, overlap.yticks),
                     labels = c(0, overlap.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())

print(overlap.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("overlap_comparison_", analysis.file.name, ".png"), overlap.figure)

###### Network analysis Hits -----

#Create data frame for graph - pVal and shared hits
analysis.type <- "HC_NA"
analysis.name <- "Network Analysis"
analysis.file.name <- "Network"


#Create levels for factors
Levels <- c("High Scoring Hits", "Post Validation Hits", analysis.name)

hits.pval <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = analysis.name
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              )),
                   data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = "High Scoring Hits"
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              ))
                   ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                               , "HitSelection" = "Post Validation Hits"
                               , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                               )))
#Introduce linebreaks
levels(hits.pval$Comparison) <- gsub(" ", "\n", levels(hits.pval$Comparison))
levels(hits.pval$Comparison) <- gsub("-", "+", levels(hits.pval$Comparison))
#Order factor of groups
hits.pval$HitSelection <- factor(hits.pval$HitSelection, levels = Levels)

#Shared hits wih pval

hits.shared <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                , "HitSelection" = analysis.name
                                , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "High Scoring Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "Post Validation Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 )))

#Introduce linebreaks
levels(hits.shared$Comparison) <- gsub(" ", "\n", levels(hits.shared$Comparison))
levels(hits.shared$Comparison) <- gsub("-", "+", levels(hits.shared$Comparison))

hits.shared$HitSelection <- factor(hits.shared$HitSelection, levels = Levels)

#Assign star designation for significance
hits.shared$star <- ""
hits.shared$star[hits.shared$Significance > .05 & hits.shared$HitSelection == analysis.name]  <- "ns"
hits.shared$star[hits.shared$Significance <= .05 & hits.shared$HitSelection == analysis.name]  <- "*"
hits.shared$star[hits.shared$Significance <= .01 & hits.shared$HitSelection == analysis.name]  <- "**"




## Create figure



pval.figure <- ggplot(hits.pval, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection) +
  ylab("-Log p Value")+
  scale_color_manual(values = pval.border)+
  scale_fill_manual(breaks = Levels,
                    values = pval.plot.colors)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  scale_y_continuous(limits=c(0, pval.max), breaks = c(0, pval.yticks),
                     labels = c(0, pval.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())
print(pval.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("pval_comparison_", analysis.file.name, ".png"), pval.figure)


overlap.figure <- ggplot(hits.shared, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection) +
  ylab("Shared Hits")+
  scale_color_manual(values = overlap.border)+
  scale_fill_manual(breaks = Levels,
                    values = overlap.plot.colors)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  geom_text(aes(label=star), colour = "black", position = position_dodge(width = 0.85), vjust=-0.15, size=10) +
  scale_y_continuous(limits=c(0, overlap.max), breaks = c(0, overlap.yticks),
                     labels = c(0, overlap.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())

print(overlap.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("overlap_comparison_", analysis.file.name, ".png"), overlap.figure)

###### Serial analysis Hits -----

#Create data frame for graph - pVal and shared hits
analysis.type <- "SERIAL_PA2t"
analysis.name <- "Serial Analysis"
analysis.file.name <- "serial"


#Create levels for factors
Levels <- c("High Scoring Hits", "Post Validation Hits", analysis.name)

hits.pval <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = analysis.name
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              )),
                   data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = "High Scoring Hits"
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              ))
                   ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                               , "HitSelection" = "Post Validation Hits"
                               , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                               )))
#Introduce linebreaks
levels(hits.pval$Comparison) <- gsub(" ", "\n", levels(hits.pval$Comparison))
levels(hits.pval$Comparison) <- gsub("-", "+", levels(hits.pval$Comparison))
#Order factor of groups
hits.pval$HitSelection <- factor(hits.pval$HitSelection, levels = Levels)

#Shared hits wih pval

hits.shared <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                , "HitSelection" = analysis.name
                                , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "High Scoring Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "Post Validation Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 )))

#Introduce linebreaks
levels(hits.shared$Comparison) <- gsub(" ", "\n", levels(hits.shared$Comparison))
levels(hits.shared$Comparison) <- gsub("-", "+", levels(hits.shared$Comparison))

hits.shared$HitSelection <- factor(hits.shared$HitSelection, levels = Levels)

#Assign star designation for significance
hits.shared$star <- ""
hits.shared$star[hits.shared$Significance > .05 & hits.shared$HitSelection == analysis.name]  <- "ns"
hits.shared$star[hits.shared$Significance <= .05 & hits.shared$HitSelection == analysis.name]  <- "*"
hits.shared$star[hits.shared$Significance <= .01 & hits.shared$HitSelection == analysis.name]  <- "**"




## Create figure



pval.figure <- ggplot(hits.pval, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection) +
  ylab("-Log p Value")+
  scale_color_manual(values = pval.border)+
  scale_fill_manual(breaks = Levels,
                    values = pval.plot.colors)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  scale_y_continuous(limits=c(0, pval.max), breaks = c(0, pval.yticks),
                     labels = c(0, pval.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())
print(pval.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("pval_comparison_", analysis.file.name, ".png"), pval.figure)


overlap.figure <- ggplot(hits.shared, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection) +
  ylab("Shared Hits")+
  scale_color_manual(values = overlap.border)+
  scale_fill_manual(breaks = Levels,
                    values = overlap.plot.colors)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  geom_text(aes(label=star), colour = "black", position = position_dodge(width = 0.85), vjust=-0.15, size=10) +
  scale_y_continuous(limits=c(0, overlap.max), breaks = c(0, overlap.yticks),
                     labels = c(0, overlap.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())

print(overlap.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("overlap_comparison_", analysis.file.name, ".png"), overlap.figure)

###### Iterative analysis Hits -----

#Create data frame for graph - pVal and shared hits
analysis.type <- "TRIAGE_PA2t"
analysis.name <- "Iterative Analysis"
analysis.file.name <- "iterative"


#Create levels for factors
Levels <- c("High Scoring Hits", "Post Validation Hits", analysis.name)

hits.pval <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = analysis.name
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              )),
                   data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = "High Scoring Hits"
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              ))
                   ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                               , "HitSelection" = "Post Validation Hits"
                               , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                               )))
#Introduce linebreaks
levels(hits.pval$Comparison) <- gsub(" ", "\n", levels(hits.pval$Comparison))
levels(hits.pval$Comparison) <- gsub("-", "+", levels(hits.pval$Comparison))
#Order factor of groups
hits.pval$HitSelection <- factor(hits.pval$HitSelection, levels = Levels)

#Shared hits wih pval

hits.shared <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                , "HitSelection" = analysis.name
                                , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "High Scoring Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "Post Validation Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 )))

#Introduce linebreaks
levels(hits.shared$Comparison) <- gsub(" ", "\n", levels(hits.shared$Comparison))
levels(hits.shared$Comparison) <- gsub("-", "+", levels(hits.shared$Comparison))

hits.shared$HitSelection <- factor(hits.shared$HitSelection, levels = Levels)

#Assign star designation for significance
hits.shared$star <- ""
hits.shared$star[hits.shared$Significance > .05 & hits.shared$HitSelection == analysis.name]  <- "ns"
hits.shared$star[hits.shared$Significance <= .05 & hits.shared$HitSelection == analysis.name]  <- "*"
hits.shared$star[hits.shared$Significance <= .01 & hits.shared$HitSelection == analysis.name]  <- "**"




## Create figure



pval.figure <- ggplot(hits.pval, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection) +
  ylab("-Log p Value")+
  scale_color_manual(values = pval.border)+
  scale_fill_manual(breaks = Levels,
                    values = pval.plot.colors)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  scale_y_continuous(limits=c(0, pval.max), breaks = c(0, pval.yticks),
                     labels = c(0, pval.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())
print(pval.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("pval_comparison_", analysis.file.name, ".png"), pval.figure)


overlap.figure <- ggplot(hits.shared, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection) +
  ylab("Shared Hits")+
  scale_color_manual(values = overlap.border)+
  scale_fill_manual(breaks = Levels,
                    values = overlap.plot.colors)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  geom_text(aes(label=star), colour = "black", position = position_dodge(width = 0.85), vjust=-0.15, size=10) +
  scale_y_continuous(limits=c(0, overlap.max), breaks = c(0, overlap.yticks),
                     labels = c(0, overlap.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())

print(overlap.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("overlap_comparison_", analysis.file.name, ".png"), overlap.figure)


###### Iterative PA two tiers analysis Hits -----

#Create data frame for graph - pVal and shared hits
analysis.type <- "TRIAGE"
analysis.name <- "Iterative Analysis"
analysis.file.name <- "iterative_noPA2t"


#Create levels for factors
Levels <- c("High Scoring Hits", "Post Validation Hits", analysis.name)

hits.pval <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = analysis.name
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              )),
                   data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = "High Scoring Hits"
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              ))
                   ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                               , "HitSelection" = "Post Validation Hits"
                               , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                               )))
#Introduce linebreaks
levels(hits.pval$Comparison) <- gsub(" ", "\n", levels(hits.pval$Comparison))
levels(hits.pval$Comparison) <- gsub("-", "+", levels(hits.pval$Comparison))
#Order factor of groups
hits.pval$HitSelection <- factor(hits.pval$HitSelection, levels = Levels)

#Shared hits wih pval

hits.shared <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                , "HitSelection" = analysis.name
                                , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "High Scoring Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "Post Validation Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 )))

#Introduce linebreaks
levels(hits.shared$Comparison) <- gsub(" ", "\n", levels(hits.shared$Comparison))
levels(hits.shared$Comparison) <- gsub("-", "+", levels(hits.shared$Comparison))

hits.shared$HitSelection <- factor(hits.shared$HitSelection, levels = Levels)

#Assign star designation for significance
hits.shared$star <- ""
hits.shared$star[hits.shared$Significance > .05 & hits.shared$HitSelection == analysis.name]  <- "ns"
hits.shared$star[hits.shared$Significance <= .05 & hits.shared$HitSelection == analysis.name]  <- "*"
hits.shared$star[hits.shared$Significance <= .01 & hits.shared$HitSelection == analysis.name]  <- "**"




## Create figure



pval.figure <- ggplot(hits.pval, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection) +
  ylab("-Log p Value")+
  scale_color_manual(values = pval.border)+
  scale_fill_manual(breaks = Levels,
                    values = pval.plot.colors)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  scale_y_continuous(limits=c(0, pval.max), breaks = c(0, pval.yticks),
                     labels = c(0, pval.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())
print(pval.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("pval_comparison_", analysis.file.name, ".png"), pval.figure)


overlap.figure <- ggplot(hits.shared, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection) +
  ylab("Shared Hits")+
  scale_color_manual(values = overlap.border)+
  scale_fill_manual(breaks = Levels,
                    values = overlap.plot.colors)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  geom_text(aes(label=star), colour = "black", position = position_dodge(width = 0.85), vjust=-0.15, size=10) +
  scale_y_continuous(limits=c(0, overlap.max), breaks = c(0, overlap.yticks),
                     labels = c(0, overlap.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())

print(overlap.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("overlap_comparison_", analysis.file.name, ".png"), overlap.figure)

#Create data frame for graph - pVal and shared hits
analysis.type <- "TRIAGE"
analysis.name <- "Iterative Analysis"
analysis.file.name <- "iterative_noPA2t"


#Create levels for factors
Levels <- c("High Scoring Hits", "Post Validation Hits", analysis.name)

hits.pval <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = analysis.name
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              )),
                   data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = "High Scoring Hits"
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              ))
                   ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                               , "HitSelection" = "Post Validation Hits"
                               , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                               )))
#Introduce linebreaks
levels(hits.pval$Comparison) <- gsub(" ", "\n", levels(hits.pval$Comparison))
levels(hits.pval$Comparison) <- gsub("-", "+", levels(hits.pval$Comparison))
#Order factor of groups
hits.pval$HitSelection <- factor(hits.pval$HitSelection, levels = Levels)

#Shared hits wih pval

hits.shared <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                , "HitSelection" = analysis.name
                                , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "High Scoring Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "Post Validation Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 )))

#Introduce linebreaks
levels(hits.shared$Comparison) <- gsub(" ", "\n", levels(hits.shared$Comparison))
levels(hits.shared$Comparison) <- gsub("-", "+", levels(hits.shared$Comparison))

hits.shared$HitSelection <- factor(hits.shared$HitSelection, levels = Levels)

#Assign star designation for significance
hits.shared$star <- ""
hits.shared$star[hits.shared$Significance > .05 & hits.shared$HitSelection == analysis.name]  <- "ns"
hits.shared$star[hits.shared$Significance <= .05 & hits.shared$HitSelection == analysis.name]  <- "*"
hits.shared$star[hits.shared$Significance <= .01 & hits.shared$HitSelection == analysis.name]  <- "**"




## Create figure



pval.figure <- ggplot(hits.pval, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection) +
  ylab("-Log p Value")+
  scale_color_manual(values = pval.border)+
  scale_fill_manual(breaks = Levels,
                    values = pval.plot.colors)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  scale_y_continuous(limits=c(0, pval.max), breaks = c(0, pval.yticks),
                     labels = c(0, pval.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())
print(pval.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("pval_comparison_", analysis.file.name, ".png"), pval.figure)


overlap.figure <- ggplot(hits.shared, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection) +
  ylab("Shared Hits")+
  scale_color_manual(values = overlap.border)+
  scale_fill_manual(breaks = Levels,
                    values = overlap.plot.colors)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  geom_text(aes(label=star), colour = "black", position = position_dodge(width = 0.85), vjust=-0.15, size=10) +
  scale_y_continuous(limits=c(0, overlap.max), breaks = c(0, overlap.yticks),
                     labels = c(0, overlap.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())

print(overlap.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("overlap_comparison_", analysis.file.name, ".png"), overlap.figure)


###### Iterative reverse analysis Hits -----

#Create data frame for graph - pVal and shared hits
analysis.type <- "REVERSE_TRIAGE_PA2t"
analysis.name <- "Reverse \n Iterative Analysis"
analysis.file.name <- "iterative_reverse"


#Create levels for factors
Levels <- c("High Scoring Hits", "Post Validation Hits", "Iterative: Pathway to Network", analysis.name)

hits.pval <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = analysis.name
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              )),
                   data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = "Iterative: Pathway to Network"
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("TRIAGE_PA2t", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("TRIAGE_PA2t", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("TRIAGE_PA2t", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              )),
                   data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = "High Scoring Hits"
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              ))
                   ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                               , "HitSelection" = "Post Validation Hits"
                               , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                               )))
#Introduce linebreaks
levels(hits.pval$Comparison) <- gsub(" ", "\n", levels(hits.pval$Comparison))
levels(hits.pval$Comparison) <- gsub("-", "+", levels(hits.pval$Comparison))
#Order factor of groups
hits.pval$HitSelection <- factor(hits.pval$HitSelection, levels = Levels)

#Shared hits wih pval

hits.shared <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                , "HitSelection" = analysis.name
                                , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "Iterative: Pathway to Network"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("TRIAGE_PA2t", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("TRIAGE_PA2t", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("TRIAGE_PA2t", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("TRIAGE_PA2t", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("TRIAGE_PA2t", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("TRIAGE_PA2t", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "High Scoring Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "Post Validation Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 )))
#Introduce linebreaks
levels(hits.shared$Comparison) <- gsub(" ", "\n", levels(hits.shared$Comparison))
levels(hits.shared$Comparison) <- gsub("-", "+", levels(hits.shared$Comparison))

hits.shared$HitSelection <- factor(hits.shared$HitSelection, levels = Levels)

#Assign star designation for significance
hits.shared$star <- ""
hits.shared$star[hits.shared$Significance > .05 & hits.shared$HitSelection == analysis.name]  <- "ns"
hits.shared$star[hits.shared$Significance <= .05 & hits.shared$HitSelection == analysis.name]  <- "*"
hits.shared$star[hits.shared$Significance <= .01 & hits.shared$HitSelection == analysis.name]  <- "**"
hits.shared$star[hits.shared$Significance > .05 & hits.shared$HitSelection == "Iterative: Pathway to Network"]  <- "ns"
hits.shared$star[hits.shared$Significance <= .05 & hits.shared$HitSelection == "Iterative: Pathway to Network"]  <- "*"
hits.shared$star[hits.shared$Significance <= .01 & hits.shared$HitSelection == "Iterative: Pathway to Network"]  <- "**"



## Create figure



pval.figure <- ggplot(hits.pval, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection, ncol = 4) +
  ylab("-Log p Value")+
  scale_color_manual(values = pval.border2)+
  scale_fill_manual(breaks = Levels,
                    values = pval.plot.colors2)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  scale_y_continuous(limits=c(0, pval.max), breaks = c(0, pval.yticks),
                     labels = c(0, pval.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())
print(pval.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("pval_comparison_", analysis.file.name, ".png"), pval.figure)


overlap.figure <- ggplot(hits.shared, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection, ncol = 4) +
  ylab("Shared Hits")+
  scale_color_manual(values = overlap.border2)+
  scale_fill_manual(breaks = Levels,
                    values = overlap.plot.colors2)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  geom_text(aes(label=star), colour = "black", position = position_dodge(width = 0.85), vjust=-0.15, size=10) +
  scale_y_continuous(limits=c(0, overlap.max), breaks = c(0, overlap.yticks),
                     labels = c(0, overlap.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())

print(overlap.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("overlap_comparison_", analysis.file.name, ".png"), overlap.figure)

#Create data frame for graph - pVal and shared hits
analysis.type <- "TRIAGE"
analysis.name <- "Iterative Analysis"
analysis.file.name <- "iterative_noPA2t"


#Create levels for factors
Levels <- c("High Scoring Hits", "Post Validation Hits", analysis.name)

hits.pval <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = analysis.name
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              )),
                   data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = "High Scoring Hits"
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              ))
                   ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                               , "HitSelection" = "Post Validation Hits"
                               , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                               )))
#Introduce linebreaks
levels(hits.pval$Comparison) <- gsub(" ", "\n", levels(hits.pval$Comparison))
levels(hits.pval$Comparison) <- gsub("-", "+", levels(hits.pval$Comparison))
#Order factor of groups
hits.pval$HitSelection <- factor(hits.pval$HitSelection, levels = Levels)

#Shared hits wih pval

hits.shared <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                , "HitSelection" = analysis.name
                                , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "High Scoring Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "Post Validation Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 )))
#Introduce linebreaks
levels(hits.shared$Comparison) <- gsub(" ", "\n", levels(hits.shared$Comparison))
levels(hits.shared$Comparison) <- gsub("-", "+", levels(hits.shared$Comparison))

hits.shared$HitSelection <- factor(hits.shared$HitSelection, levels = Levels)

#Assign star designation for significance
hits.shared$star <- ""
hits.shared$star[hits.shared$Significance > .05 & hits.shared$HitSelection == analysis.name]  <- "ns"
hits.shared$star[hits.shared$Significance <= .05 & hits.shared$HitSelection == analysis.name]  <- "*"
hits.shared$star[hits.shared$Significance <= .01 & hits.shared$HitSelection == analysis.name]  <- "**"




## Create figure



pval.figure <- ggplot(hits.pval, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection) +
  ylab("-Log p Value")+
  scale_color_manual(values = pval.border2)+
  scale_fill_manual(breaks = Levels,
                    values = pval.plot.colors)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  scale_y_continuous(limits=c(0, pval.max), breaks = c(0, pval.yticks),
                     labels = c(0, pval.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())
print(pval.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("pval_comparison_", analysis.file.name, ".png"), pval.figure)


overlap.figure <- ggplot(hits.shared, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7, size = border.size) +    
  facet_wrap(~HitSelection) +
  ylab("Shared Hits")+
  scale_color_manual(values = overlap.border2)+
  scale_fill_manual(breaks = Levels,
                    values = overlap.plot.colors)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  geom_text(aes(label=star), colour = "black", position = position_dodge(width = 0.85), vjust=-0.15, size=10) +
  scale_y_continuous(limits=c(0, overlap.max), breaks = c(0, overlap.yticks),
                     labels = c(0, overlap.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())

print(overlap.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("overlap_comparison_", analysis.file.name, ".png"), overlap.figure)



################ FIgure 3 ################


# Get Data Frame
##-- 
setwd("~/Documents/Analysis/HIV/TRIAGEHIV_May2019/Comparison/PermutationTests")
comparison.analysis.df <- read.csv(paste0("ComparisonAnalysis_1000permTest.csv"), stringsAsFactors = FALSE)

###-- Seperate Dataframes by analysis and filter by analysis you want to include in figure

analysis.included <- c("UGP", "HC", "HC_PA", "HC_PA2t", "HC_NA", "SERIAL_PA2t", "TRIAGE_PA2t")

overlap.df <- comparison.analysis.df[which(grepl("Overlap", comparison.analysis.df$X) == T), ]
overlap.df <- comparison.analysis.df[which(comparison.analysis.df$X %in% paste0(analysis.included, "_Overlap")), ]

pval.df <- comparison.analysis.df[which(grepl("negLogP", comparison.analysis.df$X) == T), ]
pval.df <- comparison.analysis.df[which(comparison.analysis.df$X %in% paste0(analysis.included, "_negLogP")), ]

perm_test.df <- comparison.analysis.df[which(grepl("Random.perc", comparison.analysis.df$X) == T), ]
perm_test.df <- comparison.analysis.df[which(comparison.analysis.df$X %in% paste0(analysis.included, "_Random.perc")), ]


#- Filter for the analysis you want to include 
##- Find max for each dataframe
overlap.max <- (max(overlap.df[, -1]) + 1)
pval.max <- ceiling(max(pval.df[, c(-1, -5)]))


###Create dataframe with percentage of max value for each cell
overlap.percmax <- overlap.df %>%
  mutate(zhou.brass = (zhou.brass / overlap.max)) %>%
  mutate(brass.konig = (brass.konig / overlap.max)) %>%
  mutate(zhou.konig = (zhou.konig / overlap.max))

pval.percmax <- pval.df %>%
  mutate(zhou.brass = (zhou.brass / pval.max)) %>%
  mutate(brass.konig = (brass.konig / pval.max)) %>%
  mutate(zhou.konig = (zhou.konig / pval.max))

###Create Dataframe with asterisk signs for perm test
perm_test.symbol <- perm_test.df %>%
  mutate(zhou.brass = ifelse(zhou.brass <= 0.01, "**", 
                             ifelse(zhou.brass <= 0.05, "*", "ns"))) %>%
  mutate(brass.konig = ifelse(brass.konig <= 0.01, "**", 
                              ifelse(brass.konig <= 0.05, "*", "ns"))) %>%
  mutate(zhou.konig = ifelse(zhou.konig <= 0.01, "**", 
                             ifelse(zhou.konig <= 0.05, "*", "ns")))


###--- Create non zero tick marks for both graphs
overlap.yticks <- c(ceiling((overlap.max/3)*1), ceiling((overlap.max/3)*2), overlap.max)
pval.yticks <- c(ceiling((pval.max/3)*1), ceiling((pval.max/3)*2), pval.max)

# ###-- Colors
# overlap.color <- "blue"
# pval.color <- "red"

overlap.ticks.color <- "#61853d"
pval.ticks.color <- "#306880"


############# Circlize -------
circos.clear()  #Clear plot 

factors = 1:14
circos.par(start.degree = 90, gap.degree = rep(c(0.2,8), 7))

circos.initialize(factors = factors, xlim = c(0, 1))

# yang yao is __ (a long segment)
add_yang_yao = function() {
  i = ((get.cell.meta.data("sector.numeric.index") + 1) / 2 )
  circos.lines(c(0, 0.01), c(0, 1), col = "black", type = "s", lwd = 1.5)
  circos.rect(0.06,0,0.3, pval.percmax[i, "zhou.brass"], col = paste0( pval.color))   ## zhou.brass
  circos.rect(0.35,0,0.65,pval.percmax[i, "brass.konig"], col = paste0(pval.color)) ## brass.konig
  circos.rect(0.70,0,1,pval.percmax[i, "zhou.konig"], col = paste0(pval.color)) ## zhou.konig
  circos.text(0, c(0, pval.yticks/pval.max), labels = c("0", pval.yticks),
              facing = "outside", niceFacing = T,
              adj = c(-0.1, c(1, 0, 0, 0)),
              cex = 0.7,
              col = pval.ticks.color)
}

# yin yao is -- (two short segments)
add_yin_yao = function() {
  i = ((get.cell.meta.data("sector.numeric.index")) / 2 )
  circos.lines(c(0.99, 1), c(0, 1), col = "black", type = "s", lwd = 1.5)
  circos.rect(0,0,0.25, overlap.percmax[i, "zhou.brass"], col = paste0(overlap.color)) ## zhou.brass
  circos.rect(0.30,0,0.60, overlap.percmax[i, "brass.konig"], col = paste0(overlap.color)) ## brass.konig
  circos.rect(0.65,0,0.94, overlap.percmax[i, "zhou.konig"], col = paste0(overlap.color)) ## zhou.konig
  circos.text(1, c(0, overlap.yticks/overlap.max), labels = c("0", overlap.yticks),
              facing = "outside", niceFacing = T,
              adj = c(1.2, c(1, 0, 0, 0)),
              cex = 0.7,
              col = overlap.ticks.color)
  circos.text(0.125, (overlap.percmax[i, "zhou.brass"] + 0.05), perm_test.symbol[i, "zhou.brass"])
  circos.text(0.45, (overlap.percmax[i, "brass.konig"] + 0.05), perm_test.symbol[i, "brass.konig"]) 
  circos.text(0.79, (overlap.percmax[i, "zhou.konig"] + 0.05), perm_test.symbol[i, "zhou.konig"]) 
}



circos.track(ylim = c(0, 1), factors = factors, bg.border = NA,
             panel.fun = function(x, y) {
               i = get.cell.meta.data("sector.numeric.index")
               if(i %in% c(1, 3, 5, 7, 9, 11, 13)) add_yang_yao() else add_yin_yao()
             }, track.height = 0.7)
################################

##### - plot for group titles
circos.clear()  #Clear plot 

analysis.titles <- c("Post Validation Hits", "High Scoring Hits", "2 Tier Pathway Analysis", "3 Tier Pathway Analysis", "3 Tier Network Analysis", "Serial Analysis", "Iterative Analysis")

factors = 1:length(analysis.titles)
circos.par(start.degree = 90, gap.degree = 8)


circos.initialize(factors = factors, xlim = c(0, 1))

circos.track(ylim = c(0, 1), factors = factors, bg.col = "white", bg.border = NA, 
             # panel.fun = function(x, y) {
             #   circos.lines(c(0, 1), c(0, 0), col = "black", type = "s", lwd = 2)
             # },
             track.height = 0.2)

circos.trackText(x = rep(0.5, length(analysis.titles)), y = rep(0.4, length(analysis.titles)),
                 labels = analysis.titles,
                 cex = 1.4, factors = factors, col = c(rep("black", 6), "red"), font = 1, facing = "bending.inside", niceFacing = T)



##############

###-----Plot for legends

circos.clear()  #Clear plot 

factors = 1:14
circos.par(start.degree = 90, gap.degree = c(5))

circos.initialize(factors = factors, xlim = c(0, 1))

circos.track(ylim = c(0, 1), factors = factors, bg.border = NA,
             panel.fun = function(x, y) {
               #circos.lines(c(0, 0.1), c(0, 1), col = "black", type = "s", lwd = 1.5)
               circos.rect(0.06,0,0.3, 1, col = "gray65")
               circos.rect(0.35,0,0.65, 0.67, col = "gray65")
               circos.rect(0.70,0,1, 0.33, col = "gray65")
             }, track.height = 0.4)



circos.lines(c(0, 1), c(0, 0.1), col = "black", type = "s", lwd = 1.5)



circos.text(for(y in 1:length(overlap.yticks))
  1, 00., labels = c("0", "2", "3", "4"),
  facing = "outside", niceFacing = F,
  adj = c(1.2, 0.5),
  cex = 0.9,
  col = "red")

circos.trackText(x = rep(c(0.6, 0.75), 7), y = rep(0.99, 14),
                 labels = rep(c("-Log P", "Shared Hits"), 7),
                 cex = 1, factors = factors, col = "black", font = 2, facing = "bending.inside")

circos.trackText(x = rep(c(0.3, 0.5), 7), y = rep(0.99, 14),
                 labels = rep(c("-Log P", "Shared Hits"), 7),
                 cex = 0.7, factors = factors, col = "gray46", font = 2, facing = "bending.inside")

analysis.titles <- c("UGP", "HC", "HC + PA", "HC* + PA", "HC* + NA", "SERIAL", "TRIAGE")


factors = 1:length(analysis.titles)
circos.trackText(x = rep(0.5, length(analysis.titles)), y = rep(0.5, length(analysis.titles)),
                 labels = analysis.titles,
                 cex = 0.4, factors = factors, col = "#000000", font = 2, facing = "downward")


# circos.track(ylim = c(0, 1), factors = factors, bg.col = "white", track.height = 0.05)
# 
# circos.trackText(x = rep(0.5, 14), y = rep(0.6, 14),
#                  labels = rep(c("-Log P", "Shared Hits"), 7),
#                  cex = 0.8, factors = factors, col = "#000000", font = 2, facing = "bending.inside")

##### - plot for group titles
circos.clear()  #Clear plot 

analysis.titles <- c("Post Validation Hits", "High Scoring Hits", "2 Tier Pathway Analysis", "3 Pathway Analysis", "3 Tier Network Analysis", "Serial Analysis", "Iterative Analysis")

factors = 1:length(analysis.titles)
circos.par(start.degree = 90, gap.degree = 12)


circos.initialize(factors = factors, xlim = c(0, 1))

circos.track(ylim = c(0, 1), factors = factors, bg.col = "black", track.height = 0.2)

circos.trackText(x = rep(0.5, length(analysis.titles)), y = rep(0.6, length(analysis.titles)),
                 labels = analysis.titles,
                 cex = 2, factors = factors, col = "white", font = 2, facing = "bending.inside")


circos.track(ylim = c(0, 1), factors = factors, bg.border = NA,
             panel.fun = function(x, y) {
               i = get.cell.meta.data("sector.numeric.index")
               
             }, track.height = 0.3)

circos.text(x + ux(2, "mm"), y + uy(2, "mm"), labels)

circos.track(ylim = c(0, 1), factors = 1:length(analysis.titles), bg.border = NA,
             panel.fun = function(x, y) {
               i = get.cell.meta.data("sector.numeric.index")
               circos.text(x + ux(2, "mm"), y + uy(2, "mm"), labels = analysis.titles[i])
             }, track.height = 0.1)


#################### Supp Figure 1 ----

triage.gen.wd <- "~/Documents/Analysis/HIV/TRIAGEHIV_May2019"

#- redefine strings from above as necessary
screen.names <- c("zhou", "brass", "konig")

PathwayDatabase <- c("KEGG"
                     #, "KEGG.hiv"
)

tiers <- c("normalized", "validated")



#- Get the TRIAGEoutput_ALL file for all screens and create data frame for iterations and hit size
for (c in 1:length(PathwayDatabase)) {
  PathwayDB <- PathwayDatabase[c]
  
  for (t in 1:length(tiers)) {
    tier <- tiers[t]
    allscreens.iterations.df <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(allscreens.iterations.df) <- c("Screen", "Iteration", "Hits")
    
    
    for (s in 1:length(screen.names)){
      screen.name <- screen.names[s]
      
      setwd(paste0(triage.gen.wd, "/TRIAGE/TRIAGEoutput/temp/"))
      
      triageoutput.all <- read.csv(paste0(screen.name, "_", tier, "_", PathwayDB, "_TRIAGEoutput_ALL.csv"), stringsAsFactors = F)
      
      #Get index of iterarion columns
      iterations.index <- c(which(grepl("ConfidenceCategory", colnames(triageoutput.all))),
                            which(grepl("Network.class.iteration", colnames(triageoutput.all))))
      
      iteration.setsize <- NULL
      for (i in 1:length(iterations.index)){
        iteration.setsize <- c(iteration.setsize, length(which(triageoutput.all[[iterations.index[i]]] == "HighConf")))
      }
      
      iterations.df <- data.frame("Screen" = screen.name, "Iteration" = 0:(length(iterations.index)-1), "Hits" = iteration.setsize)
      
      
      allscreens.iterations.df <- rbind(allscreens.iterations.df, iterations.df)
      
      #- assing specific TRIAGE output file to variable
      assign(paste0(screen.name, ".", tier, ".", PathwayDB, ".triageoutput.all"), triageoutput.all)
    }
    
    #- Create figure
    iterations.plot <-ggplot(allscreens.iterations.df, aes(x = Iteration, y = Hits, group= Screen)) +
      geom_line(aes(color=Screen), size = 1)+
      geom_point(aes(color=Screen), size = 3)+
      scale_x_continuous(name="Iterations of TRIAGE Analysis", breaks = 0:(max(allscreens.iterations.df$Iteration)+1), limits=c(0, (max(allscreens.iterations.df$Iteration)))) +
      scale_y_continuous(name="Hits Selected", limits=c((min(allscreens.iterations.df$Hits) - 50), (max(allscreens.iterations.df$Hits) + 50)))+
      theme_linedraw()+
      ggtitle(paste0("Iterations and Hit Sets by TRIAGE for ", tier," HDF Screens"))+
      theme(plot.title = element_text(lineheight=3, family="Times", color="blue", size=15, hjust = 0.5, vjust = 1))+
      theme(axis.title.y = element_text(size = rel(0.8), angle = 90, vjust = 2, face = "bold"))+ 
      theme(axis.title.x = element_text(size = rel(0.8), angle = 00, vjust = -1, face = "bold"))+
      # theme(text = element_text(size=15),
      #       axis.text.x = element_text(angle=00, hjust=0.5),
      #       panel.grid.minor.y = element_line(colour = "grey",size=0.5, linetype = "dashed"),
      #       panel.grid.minor.x = element_blank(),
      #       panel.grid.major.y = element_line(colour = "grey",size=0.5, linetype = "dashed"),
      #       panel.grid.major.x = 
      #         element_blank()
      #         #element_line(colour = "grey", size = 0.75)
      #       )
      theme(text = element_text(size=15),
            axis.text.x = element_text(angle=00, hjust=0.5),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = 
              element_line(colour = "black",size=0.5, linetype = "dotted")
            #element_line(colour = "grey", size = 0.75)
      )
    
    
    
    print(iterations.plot)
    
    setwd(paste0(results.wd, "EnrichmentFigures/"))
    
    save_plot(paste0( tier, PathwayDB, "triageiterations.png"), iterations.plot, base_height = 5, base_aspect_ratio = 2)
    
  }
}


############ Supp Figure 3 ----------

###### Iterative reverse analysis Hits -----

#Create data frame for graph - pVal and shared hits
analysis.type <- "TRIAGE_PA2t"
analysis.name <-"TRIAGE: High Scoring"
analysis.file.name <- "post_val_triage"


#Create levels for factors
Levels <- c("High Scoring Hits", "Post Validation Hits", analysis.name, "TRIAGE: Post Validation")

hits.pval <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = analysis.name
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              )),
                   data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = "TRIAGE: Post Validation"
                              , "LogpVal" = c(log10(2.8e-10)*-1
                                              ,(log10(4.1e-11)*-1)
                                              ,(log10(2.8e-11)*-1))
                   ),
                   data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                              , "HitSelection" = "High Scoring Hits"
                              , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                              ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                              ))
                   ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                               , "HitSelection" = "Post Validation Hits"
                               , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                               ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_negLogP")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                               )))
#Introduce linebreaks
levels(hits.pval$Comparison) <- gsub(" ", "\n", levels(hits.pval$Comparison))
levels(hits.pval$Comparison) <- gsub("-", "+", levels(hits.pval$Comparison))
#Order factor of groups
hits.pval$HitSelection <- factor(hits.pval$HitSelection, levels = Levels)

#Shared hits wih pval

hits.shared <- rbind(data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                , "HitSelection" = analysis.name
                                , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                     ,comparison.analysis.df[which(comparison.analysis.df$X == paste0(analysis.type, "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "TRIAGE: Post Validation"
                                 , "LogpVal" = c(34, 33, 29)
                                 , "Significance" = c(1,1,1)
                     )
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "High Scoring Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("HC", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 ))
                     ,data.frame("Comparison" = c("Zhou - Brass", "Brass - Konig", "Zhou - Konig")
                                 , "HitSelection" = "Post Validation Hits"
                                 , "LogpVal" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                 ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Overlap")), which(colnames(comparison.analysis.df) == "zhou.konig")])
                                 , "Significance" = c(comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.brass")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "brass.konig")]
                                                      ,comparison.analysis.df[which(comparison.analysis.df$X == paste0("UGP", "_Random.perc")), which(colnames(comparison.analysis.df) == "zhou.konig")]
                                 )))
#Introduce linebreaks
levels(hits.shared$Comparison) <- gsub(" ", "\n", levels(hits.shared$Comparison))
levels(hits.shared$Comparison) <- gsub("-", "+", levels(hits.shared$Comparison))

hits.shared$HitSelection <- factor(hits.shared$HitSelection, levels = Levels)

#Assign star designation for significance
hits.shared$star <- ""
hits.shared$star[hits.shared$Significance > .05 & hits.shared$HitSelection == analysis.name]  <- "ns"
hits.shared$star[hits.shared$Significance <= .05 & hits.shared$HitSelection == analysis.name]  <- "*"
hits.shared$star[hits.shared$Significance <= .01 & hits.shared$HitSelection == analysis.name]  <- "**"
hits.shared$star[hits.shared$Significance > .05 & hits.shared$HitSelection == "Iterative: Pathway to Network"]  <- "ns"
hits.shared$star[hits.shared$Significance <= .05 & hits.shared$HitSelection == "Iterative: Pathway to Network"]  <- "*"
hits.shared$star[hits.shared$Significance <= .01 & hits.shared$HitSelection == "Iterative: Pathway to Network"]  <- "**"



## Create figure



pval.figure <- ggplot(hits.pval, aes(y=LogpVal, x=Comparison, fill=HitSelection, color=HitSelection)) + 
  geom_bar( stat="identity", width = 0.7) +    
  facet_wrap(~HitSelection, ncol = 4) +
  ylab("-Log p Value")+
  scale_color_manual(values = pval.border2)+
  scale_fill_manual(breaks = Levels,
                    values = pval.plot.colors2)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  scale_y_continuous(limits=c(0, pval.max), breaks = c(0, pval.yticks),
                     labels = c(0, pval.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())
print(pval.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("pval_comparison_", analysis.file.name, ".png"), pval.figure)


overlap.figure <- ggplot(hits.shared, aes(y=LogpVal, x=Comparison, fill=HitSelection, color = HitSelection)) + 
  geom_bar( stat="identity", width = 0.7) +    
  facet_wrap(~HitSelection, ncol = 4) +
  ylab("Shared Hits")+
  scale_color_manual(values = overlap.border2)+
  scale_fill_manual(breaks = Levels,
                    values = overlap.plot.colors2)+
  #scale_color_continuous(breaks = Levels,
  #values = pval.plot.colors)+
  #geom_text(aes(label=star), colour = "black", position = position_dodge(width = 0.85), vjust=-0.15, size=10) +
  scale_y_continuous(limits=c(0, overlap.max), breaks = c(0, overlap.yticks),
                     labels = c(0, overlap.yticks)) +
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 2/1,
        panel.background = element_blank())

print(overlap.figure )
setwd(paste0(results.wd, "EnrichmentFigures/"))

ggsave(paste0("overlap_comparison_", analysis.file.name, ".png"), overlap.figure)



############# Q1: How is the RAW data related to the Normalized Data -----

#libraries
triage.gen.wd <- "~/Documents/Analysis/HIV/TRIAGEHIV_May2019"
hiv.raw.wd <- paste0(triage.gen.wd, "/OriginalData/AuthorProvided") 

# get files ---> for each screen get name.raw.df, name.bhaskar.df, name.card.df

#Zhou
zhou.raw.df <- read.csv(paste0(hiv.raw.wd, "/Zhou/MerckDataOriginal_genes.csv"), stringsAsFactors = F)

#Brass
brass.raw.df <- read.csv(paste0(hiv.raw.wd, "/Brass/SMARTpoolHIVscreenforIainFeb26.csv"), stringsAsFactors = F)

#Konig
konig.raw.df <- read.csv(paste0(hiv.raw.wd, "/Konig/ChandaOriginalData.csv"), stringsAsFactors = F)

#####- Normalization
#Zhou-
#Log transfre inhibition to expression levels and Zscore for viability and remove non human genes
#Define humna gene list
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)

zhou.normalized.df <-  zhou.raw.df %>%
  mutate(X48Hexpression = log(1 - (0.01 * X48HpercINH))) %>%
  mutate(X96Hexpression = log(1 - (0.01 * X96HpercINH))) %>%
  subset(EntrezID %in% mapped_genes) %>%
  mutate(X48hViab.Zscore = as.numeric(scale(percViability48H, center = TRUE, scale = TRUE))) %>%
  mutate(X96hViab.Zscore = as.numeric(scale(percViability96H, center = TRUE, scale = TRUE)))

#colnames
zhou.raw.name <- "X48HpercINH"
zhou.normalized.name <- "X48Hexpression"
zhou.confidence.name <- "X48hViab.Zscore"

#Brass-

# - Remove rows with no Entrez, log normalize cell count.
brass.normalized.df <- brass.raw.df %>%
  subset(is.na(EntrezID) == F) %>%
  subset(EntrezID != 0) %>%
  mutate(Log10CellNumber = log10(NormalizedCellNumber)) %>%
  mutate(UniqueID = paste0(Plate, Well)) #Unique identifier



#-nromalize percent infected and cell count by plate 
plates <- unique(brass.normalized.df$Plate)   #Get list of plate ID

brass.normalized.df$PercInfected.Zscore <- NA #Name column for normalized values

PlateNormCol <- which(colnames(brass.normalized.df) == "PercInfected.Zscore")
DataCol <- which(colnames(brass.normalized.df) == "NormalizedPercentInfected")

for(i in 1:length(plates)){
  indPlate <- which(brass.normalized.df$Plate == plates[i])
  
  brass.normalized.df[indPlate, PlateNormCol] <- (brass.normalized.df[indPlate,DataCol] -  mean(brass.normalized.df[indPlate,DataCol],na.rm = TRUE))/
    sd(brass.normalized.df[indPlate,DataCol],na.rm = TRUE)
}

brass.normalized.df$CellNumber.Zscore <- NA

PlateNormCol <- which(colnames(brass.normalized.df) == "CellNumber.Zscore")
DataCol <- which(colnames(brass.normalized.df) == "Log10CellNumber")

for(i in 1:length(plates)){
  indPlate <- which(brass.normalized.df$Plate == plates[i])
  
  brass.normalized.df[indPlate, PlateNormCol] <- (brass.normalized.df[indPlate,DataCol] -  mean(brass.normalized.df[indPlate,DataCol],na.rm = TRUE))/
    sd(brass.normalized.df[indPlate,DataCol],na.rm = TRUE)
}


brass.raw.name <- "NormalizedPercentInfected"
brass.normalized.name <- "PercInfected.Zscore"
brass.confidence.name <- "CellNumber.Zscore"

#- Konig
konig.normalized.df <- konig.raw.df %>%
  tidyr::separate(Well_ID, c("Plate", "Well") ) %>%
  mutate(PlateID = paste(Plate, substring(Library, 4), sep = "")) %>%
  mutate(negLogP = LogP*-1) %>%
  mutate(LogScore = log(Score)) %>%
  mutate(UniqueID = paste0(PlateID, "_", Well)) #Unique identifier


konig.normalized.df$Plate <- NULL  

names(konig.normalized.df)[which(colnames(konig.normalized.df) == "Gene_ID")] <- "EntrezID"

#-nromalize percent infected and cell count by plate 
plates <- unique(konig.normalized.df$PlateID)   #Get list of plate ID

konig.normalized.df$Zscore <- NA #Name column for normalized values

PlateNormCol <- which(colnames(konig.normalized.df) == "Zscore")
DataCol <- which(colnames(konig.normalized.df) == "LogScore")

for(i in 1:length(plates)){
  indPlate <- which(konig.normalized.df$PlateID == plates[i])
  
  konig.normalized.df[indPlate, PlateNormCol] <- (konig.normalized.df[indPlate,DataCol] -  mean(konig.normalized.df[indPlate,DataCol],na.rm = TRUE))/
    sd(konig.normalized.df[indPlate,DataCol],na.rm = TRUE)
}

konig.raw.name <- "Score"
konig.normalized.name <- "Zscore"
konig.confidence.name <- "negLogP"


#######Q1 Figures
screen.names <- c("zhou", "brass", "konig")
raw.columns <- c(zhou.raw.name, brass.raw.name, konig.raw.name)
normalized.columns <- c(zhou.normalized.name, brass.normalized.name, konig.normalized.name)
confidence.columns <- c(zhou.confidence.name, brass.confidence.name, konig.confidence.name)
confidence.names <- c("Viability", "Viability", "-LogP")

# Original vs. Normalized
for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  
  assign(paste0(screen.name, ".plot"), ggplot(get(paste0(screen.name, ".normalized.df")))+
           geom_point(aes(x = get(normalized.columns[s]), y = get(raw.columns[s])))+
           scale_shape_manual()+
           xlab('Zscore') + ylab('Published Scores')+
           theme(legend.position="none")+
           ggtitle(paste0(toupper(screen.name), ": Score Transformation")))
  
  setwd(paste0(triage.gen.wd, "/OriginalData/AuthorProvided/Figures"))
  save_plot(paste0(screen.name, ".scoretransform.png"), get(paste0(screen.name, ".plot")))
}

# Normalized vs. Confidence
for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  
  assign(paste0(screen.name, ".plot"), ggplot(get(paste0(screen.name, ".normalized.df")))+
           geom_point(aes(x = get(normalized.columns[s]), y = get(confidence.columns[s])))+
           scale_shape_manual()+
           xlab('Zscore') + ylab(confidence.names[s])+
           theme(legend.position="none")+
           ggtitle(paste0(toupper(screen.name), ": Score & Confidence")))
  
  setwd(paste0(triage.gen.wd, "/OriginalData/AuthorProvided/Figures"))
  save_plot(paste0(screen.name, ".scoreandconf.png"), get(paste0(screen.name, ".plot")))
}


# count versus score
for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  
  assign(paste0(screen.name, ".plot"), ggplot(get(paste0(screen.name, ".normalized.df")))+
           geom_histogram(aes(x = get(normalized.columns[s])), binwidth = 0.1, color="black", fill="black")+
           ylab("Frequency") + xlab("Zscore")+
           theme(legend.position="none")+
           ggtitle(paste0(toupper(screen.name), ": Score Distribution")))
  
  setwd(paste0(triage.gen.wd, "/OriginalData/AuthorProvided/Figures"))
  save_plot(paste0(screen.name, ".scorefrequency.png"), get(paste0(screen.name, ".plot")))
  
}

############# Q2: Distribution of median scores for each study -----
#- get median zscore and median confidence measure for each screen
for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  #aggregate the score
  agg.score <- aggregate(get(normalized.columns[s]) ~ EntrezID, get(paste0(screen.name, ".normalized.df")), median)
  names(agg.score)[2] <- paste0((normalized.columns[s]))
  agg.conf <- aggregate(get(confidence.columns[s]) ~ EntrezID, get(paste0(screen.name, ".normalized.df")), median)
  names(agg.conf)[2] <- paste0((confidence.columns[s]))
  
  assign(paste0(screen.name, ".median.df"), merge(agg.score, agg.conf))
  
  
}

#Q2 figures
for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  
  assign(paste0(screen.name, ".plot"), ggplot(get(paste0(screen.name, ".median.df")))+
           geom_point(aes(x = get(normalized.columns[s]), y = get(confidence.columns[s])))+
           scale_shape_manual()+
           xlab('Median Zscore') + ylab(paste0("Median",confidence.names[s]))+
           theme(legend.position="none")+
           ggtitle(paste0(toupper(screen.name), ": Median Score & Confidence")))
  
  setwd(paste0(triage.gen.wd, "/OriginalData/AuthorProvided/Figures"))
  save_plot(paste0(screen.name, ".median.png"), get(paste0(screen.name, ".plot")))
}

############# Q3: Where are the hits selected for follow up in each study -----
# #Zhou - 
# #- to find the selected hits use the data frame and you can see which row has values in the "Confirm1" column indicating that it was followed up on,
# ##!!!! Use Compound instead of ENtrezID as that is what it is unique
# zhou.followup.matrix <- zhou.normalized.df$Compound[-c(which(is.na(zhou.normalized.df$X48HpercINH_Confirm1)), which(is.na(zhou.normalized.df$EntrezID)))]
# zhou.uniqueID <- "Compound"
# 
# #Brass - Using the criteria mentioned in the paper, with caveat that we only have data for p24 measure, not for beta gal measurement
# brass.followup.matrix <- brass.normalized.df$UniqueID[intersect(which(brass.normalized.df$PercInfected.Zscore <= -2), which(brass.normalized.df$CellNumber.Zscore >= -2))] 
# brass.uniqueID <- "UniqueID"
# 
# #Konig - 
# #Create 2 siRNA. above cutoff
# zscore.threshold <- -1
# 
# konig.categories2in1 <- group_by(konig.normalized.df, EntrezID) %>%
#   summarise(TwosiRNA.1 = 
#               ifelse(sort(Zscore, decreasing = F)[2] <= zscore.threshold   #2 gRNAs above threshold
#                      , 1, 0))
# #### Now add the Entrez ID specific log P
# negLogP.agg <- aggregate(negLogP ~ EntrezID, konig.normalized.df, min)
# 
# #Combine with Two siRNA below -1 data frame
# Agg.merge <- merge(konig.categories2in1, negLogP.agg)


############# Q4: Where are the confirmed hits in each study -----
#Get list of hits for each study
Supplement.WD <- "/Users/sakatz/Documents/Analysis/HIV/SupplemntaryTables/"

###  Get Author Selected Hits and create matrices 
# Number of hits per screen
zhou.hits.count <- 232
brass.hits.count <- 281
konig.hits.count <- 295

setwd(Supplement.WD)
AuthorSelectedHits.df <- read.csv("AuthorSelectedHits_HIV.csv", stringsAsFactors = F)

zhou.top.hits <- matrix(AuthorSelectedHits.df$ZhouHits[1:zhou.hits.count])
brass.top.hits <- matrix(AuthorSelectedHits.df$BrassHits[1:brass.hits.count])
konig.top.hits <- matrix(AuthorSelectedHits.df$KonigHits[1:konig.hits.count])

#Create figures highlighting selected hits

# Normalized vs. Confidence
for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  
  assign(paste0(screen.name, ".plot"), ggplot(get(paste0(screen.name, ".median.df")))+
           geom_point(aes(x = get(normalized.columns[s]), y = get(confidence.columns[s]), colour = EntrezID %in% get(paste0(screen.name, ".top.hits"))))+
           geom_point(data = subset(get(paste0(screen.name, ".median.df")), EntrezID %in% get(paste0(screen.name, ".top.hits"))),
                      aes(x = get(normalized.columns[s]), y =get(confidence.columns[s]), colour = EntrezID %in% get(paste0(screen.name, ".top.hits"))))+
           scale_colour_manual(name = 'Selected Hits', values = setNames(c('blue','black'),c(T, F))) +
           scale_shape_manual()+
           xlab('Zscore') + ylab(confidence.names[s])+
           theme(legend.position="none")+
           ggtitle(paste0(toupper(screen.name), ": Median Score & Confidence")))
  
  setwd(paste0(triage.gen.wd, "/OriginalData/AuthorProvided/Figures"))
  save_plot(paste0(screen.name, ".postvalhits.png"), get(paste0(screen.name, ".plot")))
}


############# Q5: What are the high scoring genes in each study -----
screen.confidence.thresholds <- c(-2.00, -2.00, 1.3)
#- Defines sizes of subsets
HighConf.size <- 400
MedConf.size <- 1000

#- get cutoff for 
for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  
  #- subset dataframe 
  #-- Ranking based on CARD and confidence measure
  subset.df <- filter(get(paste0(screen.name, ".median.df")), get(confidence.columns[s]) > screen.confidence.thresholds[s]) %>%
    arrange(get(normalized.columns[s]))
  
  assign(paste0(screen.name, ".scorecutoff"), subset.df[[normalized.columns[s]]][HighConf.size])
}

# Visualize

# Normalized vs. Confidence  with area highlight
for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  
  assign(paste0(screen.name, ".plot"), ggplot(get(paste0(screen.name, ".median.df")))+
           geom_rect(aes(xmin=-Inf,xmax=get(paste0(screen.name, ".scorecutoff")),ymin=screen.confidence.thresholds[s],ymax=Inf),alpha=1,fill="firebrick1")+
           geom_point(aes(x = get(normalized.columns[s]), y = get(confidence.columns[s])))+
           scale_shape_manual()+
           xlab('Zscore') + ylab(paste0(confidence.names[s]))+
           theme(legend.position="none")+
           ggtitle(paste0(toupper(screen.name), ": Median Score & Confidence")))
  
  setwd(paste0(triage.gen.wd, "/OriginalData/AuthorProvided/Figures"))
  save_plot(paste0(screen.name, ".tophitshighlight.png"), get(paste0(screen.name, ".plot")))
}

#with top hits highlight
for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  
  assign(paste0(screen.name, ".plot"), ggplot(get(paste0(screen.name, ".median.df")))+
           geom_rect(aes(xmin=-Inf,xmax=get(paste0(screen.name, ".scorecutoff")),ymin=screen.confidence.thresholds[s],ymax=Inf),alpha=1,fill="firebrick1")+
           geom_point(aes(x = get(normalized.columns[s]), y = get(confidence.columns[s]), colour = EntrezID %in% get(paste0(screen.name, ".top.hits"))))+
           geom_point(data = subset(get(paste0(screen.name, ".median.df")), EntrezID %in% get(paste0(screen.name, ".top.hits"))),
                      aes(x = get(normalized.columns[s]), y =get(confidence.columns[s]), colour = EntrezID %in% get(paste0(screen.name, ".top.hits"))))+
           scale_colour_manual(name = 'Selected Hits', values = setNames(c('blue','black'),c(T, F))) +
           scale_shape_manual()+
           xlab('Zscore') + ylab(confidence.names[s])+
           theme(legend.position="none")+
           ggtitle(paste0(toupper(screen.name), ": Median Score & Confidence")))
  
  setwd(paste0(triage.gen.wd, "/OriginalData/AuthorProvided/Figures"))
  save_plot(paste0(screen.name, ".postvalplustop.png"), get(paste0(screen.name, ".plot")))
}


############# Q6: How to segment the data for TRIAGE -----
#creat tiers

for (s in 1:length(screen.names)){
  screen.name <- screen.names[s]
  
  #- subset dataframe 
  #-- Ranking based on CARD and confidence measure
  subset.df <- filter(get(paste0(screen.name, ".median.df")), get(confidence.columns[s]) > screen.confidence.thresholds[s]) %>%
    arrange(get(normalized.columns[s])) %>% mutate(RANK.hits = row_number())
  
  #-- Ranking based on CARD and confidence measure, excluding the ones selected by authors
  subset.nohits.df <- filter(get(paste0(screen.name, ".median.df")), (get(confidence.columns[s]) > screen.confidence.thresholds[s]) & 
                               !(EntrezID %in% get(paste0(screen.name, ".top.hits")))) %>%
    arrange(get(normalized.columns[s])) %>% mutate(RANK.nonhits = row_number())
  
  #-merge to main dataframe
  assign(paste0(screen.name, ".median.df"), merge(get(paste0(screen.name, ".median.df")),
                                                  subset.df[, c("EntrezID", "RANK.hits")],
                                                  all.x = T))
  
  assign(paste0(screen.name, ".median.df"), merge(get(paste0(screen.name, ".median.df")),
                                                  subset.nohits.df[, c("EntrezID", "RANK.nonhits")],
                                                  all.x = T))
  #- remove NAs
  
  
  triaginput.df <- get(paste0(screen.name, ".median.df")) %>%
    mutate(normalized.triage = if_else(RANK.hits <= HighConf.size, 1, if_else(RANK.hits <= (HighConf.size + MedConf.size), 0.5, 0))) %>% mutate(normalized.triage = replace(normalized.triage, is.na(normalized.triage), 0)) %>%
    mutate(validated.triage = if_else(EntrezID %in%  get(paste0(screen.name, ".top.hits")), 1, if_else(RANK.nonhits <= MedConf.size, 0.5, 0))) %>% mutate(validated.triage = replace(validated.triage, is.na(validated.triage), 0))
  
  
  assign(paste0(screen.name, ".triageinput.df"), triaginput.df)
}


#-- visualize
tiers <- c("validated")

for (t in 1:length(tiers)){
  tier <- tiers[t]
  
  for (s in 1:length(screen.names)){
    screen.name <- screen.names[s]
    
    HiConf.hits <- as.matrix(get(paste0(screen.name, ".triageinput.df"))[which(get(paste0(screen.name, ".triageinput.df"))[[paste0(tier, ".triage")]] == 1), c("EntrezID")])
    MedConf.hits <- as.matrix(get(paste0(screen.name, ".triageinput.df"))[which(get(paste0(screen.name, ".triageinput.df"))[[paste0(tier, ".triage")]] == 0.5), c("EntrezID")])
    
    
    assign(paste0(screen.name, ".plot"), ggplot(get(paste0(screen.name, ".triageinput.df")))+
             geom_point(aes(x = get(normalized.columns[s]), y = get(confidence.columns[s]), colour = as.character(get(paste0(tier, ".triage")))), size=0.5)+
             geom_point(data = subset(get(paste0(screen.name, ".triageinput.df")), get(paste0(tier, ".triage")) == 0.5),
                        aes(x = get(normalized.columns[s]), y = get(confidence.columns[s]), colour = as.character(get(paste0(tier, ".triage")))), size=0.5)+
             geom_point(data = subset(get(paste0(screen.name, ".triageinput.df")), get(paste0(tier, ".triage")) == 1),
                        aes(x = get(normalized.columns[s]), y = get(confidence.columns[s]), colour = as.character(get(paste0(tier, ".triage")))), size=0.5)+
             #scale_colour_manual(name = 'Selected Hits', values = setNames(c('red', 'blue','black')), c("High Conf", "MidConf", "LowConf")) +
             scale_colour_manual(name = 'Selected Hits', values = c("1" = "red", "0.5" = "blue", "0" = "black")) +
             scale_shape_manual()+
             xlab('Zscore') + ylab(confidence.names[s])+
             theme(legend.position="none")+
             ggtitle(paste0(toupper(screen.name), ": post-", tier, " TRIAGE input")))
    
    setwd(paste0(triage.gen.wd, "/TRIAGE/TRIAGEinput/Figures"))
    save_plot(paste0(screen.name, ".", tier, "triageinput.png"), get(paste0(screen.name, ".plot")))
  }
}


