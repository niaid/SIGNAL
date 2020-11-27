#################################----
## TRIAGE Example CRISPR Screen
## and RNAseq
##################################----


# Samuel Katz
# Sept. 23, 2020

#########Load libraries ------

library("dplyr")
library("ggplot2")
library('ggthemes')
library('org.Mm.eg.db')
library('org.Hs.eg.db')
library("igraph")
library("tidyr")
library("VennDiagram")
library("dplyr")
library("igraph")
library("ggraph")
library('RColorBrewer')
library('tidyverse')
# library('ROCR')
library('caret')
#library('plotROC')
library('wesanderson')


#########Set working directories ----
parnas.wd <- "~/OneDrive\ -\ National\ Institutes\ of\ Health/TRIAGE/CRISPRscreen/Parnas/"
analysis.wd <- "~/OneDrive\ -\ National\ Institutes\ of\ Health/Analyis/"
triage.wd <- paste0(analysis.wd, "TRIAGEscripts/")
TRIAGEoutput.dir <- paste0(parnas.wd, "TRIAGEoutput/" )
triage.path <- "/Users/sakatz/OneDrive\ -\ National\ Institutes\ of\ Health/Analyis/TRIAGEscripts/"
function.wd <- "~/OneDrive\ -\ National\ Institutes\ of\ Health/TRIAGE/Functions/"
database.wd <- "~/OneDrive\ -\ National\ Institutes\ of\ Health/Databases/"



######## Import Parnas Data ---
setwd(parnas.wd)

parnas.df <- read.csv("Parnas_suppData.csv", stringsAsFactors = F)


#rename column1 GeneSymbol
# colnames(parnas.df)[which(names(parnas.df) == "Column1")] <- "GeneSymbol"


#Convert columns that are character to numeric
cols.num <- which(sapply(parnas.df, class) == "character")[-1]

parnas.df[cols.num] <- sapply(parnas.df[cols.num], as.numeric)



## Add FDR
parnas.df$logFC_fdr <- p.adjust(parnas.df$pvalue, method="hochberg")

####### Plots Parnas Data -------

#Positive Regulator Rank versus Log2FC

posRANKvFC.plot <- 
  ggplot(parnas.df, aes(y=log2FoldChange, x=Rank_positive_regulators)) + 
  geom_point()+
  ggtitle("Log FC vs. Pos Reg Rank")+
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title=element_text(size=25,face="bold"),
        axis.line.y = element_line(colour = 'black', size = 0.5),
        axis.line.x = element_line(colour = 'black', size = 0.5),
        panel.background = element_blank(),   
        #panel.grid.major = element_line(colour = "grey"), 
        #panel.grid.minor = element_line(colour = "grey", size = 0.25),
        plot.title = element_text(size = 30, hjust = 0.5, vjust = 2, 
                                  margin = margin(b = 6, t = 6)))


setwd(paste0(parnas.wd, "Figures/"))
ggsave("posRANKvFC_parnas.png", width = 12.5, height = 12.5)

# Log 2 fold change versus log2 pvalue






FCvPval.plot <- 
  ggplot(parnas.df, aes(y=-1*log2(pvalue), x=log2FoldChange)) + 
  geom_point()+
  ggtitle("Log FC vs. p Value")+
  #geom_hline(yintercept=-1*log2(0.1), linetype="solid", color = "red")+
  #geom_vline(yintercept=-1*log2(0.05), linetype="solid", color = "red")+
  #geom_hline(yintercept = -1*log2(0.05),  color = "red")+
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title=element_text(size=25,face="bold"),
        axis.line.y = element_line(colour = 'black', size = 0.5),
        axis.line.x = element_line(colour = 'black', size = 0.5),
        panel.background = element_blank(),   
        #panel.grid.major = element_line(colour = "grey"), 
        #panel.grid.minor = element_line(colour = "grey", size = 0.25),
        plot.title = element_text(size = 30, hjust = 0.5, vjust = 2, 
                                  margin = margin(b = 6, t = 6)))


setwd(paste0(parnas.wd, "Figures/"))
ggsave("FCvPval_parnas.png", width = 12.5, height = 12.5)





##################---
# Parnas Secondary Screen ----
#################---

setwd(parnas.wd)

parnas.secondary.df <- read.csv("Parnas_secondaryScreen.csv", stringsAsFactors = F) %>%
  filter(zscore != "na")

#Convert columns that are character to numeric
cols.num <- which(sapply(parnas.secondary.df, class) == "character")[-1]

parnas.secondary.df[cols.num] <- sapply(parnas.secondary.df[cols.num], as.numeric)

parnas.secondary.all.matrix <- parnas.secondary.df$GeneSymbol


#Visualize which ones were included in secondary

#Add membership to parnas.df

parnas.df <- parnas.df %>%
  mutate(in_secondary = ifelse(GeneSymbol %in% parnas.secondary.all.matrix, "1", "0"))


FCvPval.secondary.plot <- 
  ggplot(parnas.df) + 
  geom_point(aes(y=-1*log2(pvalue), x=log2FoldChange, colour = in_secondary))+
  geom_point(data = subset(parnas.df, in_secondary == "1"),  aes(y=-1*log2(pvalue), x=log2FoldChange, colour = in_secondary))+
  scale_color_manual(values = c("black", "red"))+
  #scale_fill_manual(values = c("black", "red"))+
  ggtitle("Log FC vs. p Value")+
  #geom_hline(yintercept=-1*log2(0.1), linetype="solid", color = "red")+
  #geom_vline(yintercept=-1*log2(0.05), linetype="solid", color = "red")+
  #geom_hline(yintercept = -1*log2(0.05),  color = "red")+
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title=element_text(size=25,face="bold"),
        axis.line.y = element_line(colour = 'black', size = 0.5),
        axis.line.x = element_line(colour = 'black', size = 0.5),
        panel.background = element_blank(),   
        #panel.grid.major = element_line(colour = "grey"), 
        #panel.grid.minor = element_line(colour = "grey", size = 0.25),
        plot.title = element_text(size = 30, hjust = 0.5, vjust = 2, 
                                  margin = margin(b = 6, t = 6)))


setwd(paste0(parnas.wd, "Figures/"))
ggsave("FCvPval_secondary_parnas.png", width = 12.5, height = 12.5)



# Visualize secondary screen scores



secondary.screen.plot <- 
  ggplot(parnas.secondary.df) + 
  geom_point(aes(y=-1*log2(FDR_calculated), x=zscore))+
  ggtitle("Zscore vs. p Value \n Secondary Screen")+
  #geom_hline(yintercept=-1*log2(0.1), linetype="solid", color = "red")+
  #geom_vline(yintercept=-1*log2(0.05), linetype="solid", color = "red")+
  #geom_hline(yintercept = -1*log2(0.05),  color = "red")+
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title=element_text(size=25,face="bold"),
        axis.line.y = element_line(colour = 'black', size = 0.5),
        axis.line.x = element_line(colour = 'black', size = 0.5),
        panel.background = element_blank(),   
        #panel.grid.major = element_line(colour = "grey"), 
        #panel.grid.minor = element_line(colour = "grey", size = 0.25),
        plot.title = element_text(size = 30, hjust = 0.5, vjust = 2, 
                                  margin = margin(b = 6, t = 6)))


setwd(paste0(parnas.wd, "Figures/"))
ggsave("ZSvPval_secondary_parnas.png", width = 12.5, height = 12.5)



# S curve of secondary screen
parnas.secondary.df <- parnas.secondary.df[with(parnas.secondary.df, order(zscore)), ]  #order dataframe in increasing value of zscore




parnas.secondary.df$Rank <- NULL
parnas.secondary.df$Rank <- 1:nrow(parnas.secondary.df) 

cutoff <- 0.7



s_curve.plot <- ggplot(parnas.secondary.df, aes(x=Rank, y=zscore))
s_curve.plot +
  geom_point(shape=16) +
  geom_hline(yintercept=cutoff, linetype="solid", 
              color = "blue", size=0.5)+
  # scale_y_continuous(breaks=c(min(BROAD.tlr.plus.entrez.ordered$ZHL),cutoff,0,max(BROAD.tlr.plus.entrez.ordered$ZHL)),
  #                    labels =round(c(min(BROAD.tlr.plus.entrez.ordered$ZHL),cutoff,0,max(BROAD.tlr.plus.entrez.ordered$ZHL)),1) ) +
  # scale_x_continuous(breaks = c(1, nrow(BROAD.tlr.plus.entrez.ordered) ))+
  labs(title = "Parnas et al", y = "Zscore")+
  # annotate(geom="text", x=0.1*nrow(BROAD.tlr.plus.entrez.ordered), y=0.75*cutoff, label="2.5% cutoff",
  #          color="black", size =7)+
  # annotate(geom="text", x=0.1*nrow(BROAD.tlr.plus.entrez.ordered), y=1.25*cutoff, 
  #          label=paste0(length(which(BROAD.tlr.plus.entrez.ordered$ZHL >= cutoff)), " hits"),
  #          color="black", size =7)+
  theme(axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=10),
        axis.title=element_text(size=20,face="bold"),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        panel.background = element_blank())

setwd(paste0(parnas.wd, "Figures/"))
ggsave("Scurve_secondary_parnas.png", width = 12.5, height = 12.5)



#### Add validation column

parnas.secondary.df <- parnas.secondary.df %>%
  mutate(validated1.5 = (ifelse(zscore >= 1.5, 1, 0))) %>%
  mutate(validated1.0 = (ifelse(zscore >= 1.0, 1, 0))) %>%
  mutate(validated0.75 = (ifelse(zscore >= 0.75, 1, 0))) 

val1_5.matrix <- parnas.secondary.df$GeneSymbol[which(parnas.secondary.df$validated1.5 == 1)]
val1_0.matrix <- parnas.secondary.df$GeneSymbol[which(parnas.secondary.df$validated1.0 == 1)]
val0_75.matrix <- parnas.secondary.df$GeneSymbol[which(parnas.secondary.df$validated0.75 == 1)]




#############--
# Create TRIAGE input ----
############--


# Get secondary screen candidates with primary screen values

parnas.secondary.pos.df <- parnas.df %>%
  filter(in_secondary == "1" & log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) 


#Find top hit score and lower score
percent_tophits <- 25

tophits.FC.cutoff <-  parnas.secondary.pos.df$log2FoldChange[round(percent_tophits * 0.01 * nrow(parnas.secondary.pos.df), 0)]

#skip last three as they seem to be outliers

medhits.FC.cutoff <- parnas.secondary.pos.df$log2FoldChange[nrow(parnas.secondary.pos.df) - 3]


pval.cutoff <- sort(parnas.secondary.pos.df$pvalue)[nrow(parnas.secondary.pos.df) - 3]








#### Top percentiles highlighted



# percent.cutoff.high <- 2.5
# percent.cutoff.medium <- 5
# 
# screen.dataframe.pval <- parnas.df %>%
#   filter(pvalue <= 0.05)
# 
# 
# FC.cutoff <- min(matrix(screen.dataframe.pval$log2FoldChange[which(screen.dataframe.pval$log2FoldChange > 0)]))


FCvPval.cutoff.plot <- 
  ggplot(parnas.df, aes(y=-1*log2(pvalue), x=log2FoldChange)) + 
  geom_point()+
  ggtitle("Log FC vs. p Value")+
  #geom_hline(yintercept=-1*log2(0.1), linetype="solid", color = "red")+
  geom_vline(xintercept= medhits.FC.cutoff, linetype="solid", color = "red")+
  geom_vline(xintercept= tophits.FC.cutoff, linetype="solid", color = "red")+
  geom_hline(yintercept = -1*log2(pval.cutoff),  color = "blue")+
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.title=element_text(size=25,face="bold"),
        axis.line.y = element_line(colour = 'black', size = 0.5),
        axis.line.x = element_line(colour = 'black', size = 0.5),
        panel.background = element_blank(),   
        #panel.grid.major = element_line(colour = "grey"), 
        #panel.grid.minor = element_line(colour = "grey", size = 0.25),
        plot.title = element_text(size = 30, hjust = 0.5, vjust = 2, 
                                  margin = margin(b = 6, t = 6)))


setwd(paste0(parnas.wd, "Figures/"))
ggsave("FCvPval_cutoffs_parnas.png", width = 12.5, height = 12.5)


#Create TRIAGE input

secondary.genesymbol <- parnas.secondary.pos.df$GeneSymbol

parnas.df.triageinput <- parnas.df %>%
  mutate(TRIAGEscore = ifelse(pvalue <= pval.cutoff & log2FoldChange >= tophits.FC.cutoff & GeneSymbol %in% secondary.genesymbol, 1, 
                               ifelse(GeneSymbol %in% secondary.genesymbol, 0.5, 0)))


######### TRIAGE Input -pub figures --------
ymax.medhits <- -1*log2(min(parnas.df.triageinput[which(parnas.df.triageinput$TRIAGEscore == 0.5), ][["pvalue"]]))
ymax.tophits <- -1*log2(min(parnas.df.triageinput[which(parnas.df.triageinput$TRIAGEscore == 1), ][["pvalue"]]))

xmax.tophits <- max(parnas.df.triageinput[which(parnas.df.triageinput$TRIAGEscore == 1), ][["log2FoldChange"]])
xmax.medhits <- max(parnas.df.triageinput[which(parnas.df.triageinput$TRIAGEscore == 0.5), ][["log2FoldChange"]])

x.axis.spacer <- 0.1
y.axis.spacer <- 0.5

y.breaks <- c(0, -1*log2(pval.cutoff), 10, 20, 30, 40)
y.labels <- as.character(round(y.breaks, 1))

x.breaks <- c(-2, 0, medhits.FC.cutoff, tophits.FC.cutoff, 2, 4)
x.labels <- as.character(round(x.breaks, 1))

tophits.color <- "#DF8731"
medhits.color <- "#6C44B2"

tophits.color.light <- "#efc297"
medhits.color.light <- "#ad95d7"

line.size <- 1.5
line.type <- "dashed"

FCvPval.cutoff.plot <- 
  ggplot(parnas.df.triageinput) + 
  geom_rect(aes(xmin=tophits.FC.cutoff - x.axis.spacer, xmax= xmax.tophits + x.axis.spacer, ymin=(-1*log2(pval.cutoff)), ymax=ymax.tophits + y.axis.spacer), fill = tophits.color.light, color=tophits.color.light, alpha=1, size = 0)+
  geom_rect(aes(xmin=medhits.FC.cutoff, xmax= xmax.medhits, ymin=(-1*log2(pval.cutoff)), ymax=ymax.tophits + y.axis.spacer), fill = medhits.color.light, color=medhits.color.light, alpha=1, size = 0)+
  geom_point(aes(y=-1* log2(pvalue), x=log2FoldChange, colour = in_secondary), size = 4)+
  geom_point(data = subset(parnas.df, in_secondary == "1"),  aes(y=-1*log2(pvalue), x=log2FoldChange, colour = in_secondary), size = 4)+
  scale_color_manual(values = c("grey", "black"))+
  ggtitle("Log FC vs. p Value")+
  geom_hline(yintercept=-1*log2(pval.cutoff), color = "red", linetype = line.type, size = line.size)+
  geom_segment(aes(x = medhits.FC.cutoff, y = -Inf, xend = medhits.FC.cutoff, yend = ymax.tophits + y.axis.spacer), color = "red", linetype = line.type, size = line.size)+
  geom_segment(aes(x = tophits.FC.cutoff, y = -Inf, xend = tophits.FC.cutoff, yend = ymax.tophits + y.axis.spacer), color = "red", linetype = line.type, size = line.size)+
  #geom_vline(xintercept= medhits.FC.cutoff, y linetype="solid", color = "red")+
  #geom_vline(xintercept= tophits.FC.cutoff, linetype="solid", color = "red")+
  #geom_hline(yintercept = -1*log2(pval.cutoff),  color = "blue")+
  #geom_rect(aes(xmin=tophits.FC.cutoff , xmax= xmax.tophits + x.axis.spacer, ymin=(-1*log2(pval.cutoff)) - y.axis.spacer, ymax=ymax.tophits + y.axis.spacer), fill = "white", color=tophits.color, alpha=0, size = 4.5)+
  #geom_rect(aes(xmin=medhits.FC.cutoff - x.axis.spacer, xmax= xmax.medhits, ymin=(-1*log2(pval.cutoff)) - y.axis.spacer, ymax=ymax.medhits + y.axis.spacer), fill = "white", color=medhits.color, alpha=0, size = 4.5)+
  scale_y_continuous(breaks = y.breaks,
                     labels = y.labels)+
  scale_x_continuous(breaks = x.breaks,
                     labels = x.labels)+
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -1),
        axis.line.y = element_line(colour = 'black', size = 2, lineend = "round"),
        axis.line.x = element_line(colour = 'black', size = 2, lineend = "round"),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 1/1,
        panel.background = element_blank())

print(FCvPval.cutoff.plot)

setwd(paste0(parnas.wd, "Figures/"))
ggsave("FCvPval_cutoffs_parnas_manuscript.png", width = 12.5, height = 12.5)


###################################---
######### Set TRIAGE analysis -----
###################################--

data <- parnas.df.triageinput

##### Get pathways and networks
#Pathway
setwd(triage.wd)
pathwayData <- read.csv(paste0("resources/KEGG2019_", "mouse", "_BiologicalProcesses.csv"), stringsAsFactors = F)

## Get GeneSymbols for pathway files so that they are uniform with the input file

pathwayData$GeneSymbol <- NULL

library('org.Mm.eg.db') 
x <- org.Mm.egSYMBOL


mapped_genes <- mappedkeys(x)
overlappingGenes <- intersect(mapped_genes, pathwayData$EntrezID)

xx <- as.list(x[!is.na(overlappingGenes)])
y <- unlist(xx)
y <- data.frame(GeneSymbol = y, EntrezID = names(y), row.names = NULL, stringsAsFactors=FALSE)

tempData <- merge(x=pathwayData, y=y, by="EntrezID", all.x =T)
pathwayData <- data.frame(tempData, stringsAsFactors = F)

#Get Network
string.version <- "10.5"

load(paste0(triage.path, "resources/String", string.version, ".", "mouse", ".experimental.highConf.igraph.Rdata"))
load(paste0(triage.path, "resources/String", string.version, ".", "mouse", ".database.highConf.igraph.Rdata"))
load(paste0(triage.path, "resources/String", string.version, ".", "mouse", ".experimental.midConf.igraph.Rdata"))
load(paste0(triage.path, "resources/String", string.version, ".", "mouse", ".database.midConf.igraph.Rdata"))

G <- graph.union(get(paste0("String", ".", "mouse", ".experimental.highConf.igraph"))
                 , get(paste0("String", ".", "mouse", ".database.highConf.igraph"))
                 , get(paste0("String", ".", "mouse", ".experimental.midConf.igraph"))
                 , get(paste0("String", ".", "mouse", ".database.midConf.igraph"))
                 
)




#Selected_STRINGnetwork.igraph <- G
message("Networks Loaded")

#message(networkType)
use.only.commnected.components <- c('Yes')

## Name of analysis for saved files
AnalysisTitle <- "parnas_crispr"

print(AnalysisTitle)


####################
# Define Data Tiers
####################

data <- data %>%
  mutate(ConfidenceCategory = ifelse(TRIAGEscore == 1, "HighConf", ifelse(TRIAGEscore == 0.5, "MedConf", "LowConf")))


### Add EntrezID





library('org.Mm.eg.db')
x <- org.Mm.egSYMBOL2EG


mapped_genes <- mappedkeys(x)
overlappingGenes <- intersect(as.character(as.list(mapped_genes)), as.character(data$GeneSymbol))

xx <- as.list(x[overlappingGenes])
y <- unlist(xx)
y <- data.frame(GeneSymbol = names(y), EntrezID = y, row.names = NULL, stringsAsFactors=FALSE)

tempData <- merge(x=data, y=y, by ="GeneSymbol")

data <- data.frame(tempData, stringsAsFactors = F)







##########################
##### Pre TRIAGE Analysis
##########################
#Get hits as data frames
ZscoreHits <- data[data$ConfidenceCategory == "HighConf", c("EntrezID", "GeneSymbol")]
ZscoreBackground <- data[data$ConfidenceCategory != "HighConf", c("EntrezID", "GeneSymbol")]

#Compute Enrichment

source(paste0(triage.wd, "PathwayAnalysis_Function.R"))
# pathwayData <- KEGG2017_Human_BP # pathway data is defined b species
Hits <- ZscoreHits$EntrezID
nonHits <- ZscoreBackground$EntrezID

file.name <- paste0(AnalysisTitle, "_ZscorePathwaysAnalysis.csv")
ZscorePath <- paste0(parnas.wd, "TRIAGEoutput/")
Zscore_PathwayAnalyis <- PathwayAnalysis(pathwayData, Hits, nonHits, file.name, ZscorePath)

########################
#### TRIAGE Analysis
########################
## Name of analysis for saved files


print(AnalysisTitle)

#Set temp working directory for generating enrichment files
tempWD <- paste0(parnas.wd, "TRIAGEoutput/", "temp/")
setwd(tempWD)

proxyScore <- "ConfidenceCategory"
iteration <- 1
counter <- TRUE

# Get a copy of the original list of high-confidence genes
originalHits <- data$GeneSymbol[data$ConfidenceCategory == "HighConf"]

source(paste0(triage.wd, "pathway_iteration_longform.R"), local = TRUE)

#Save original file format
data.original <- data.frame(data, stringsAsFactors = F)

# Perform iterative TRIAGE analysis
while (counter == TRUE) {
  
  ## Show progress
  
  Hits <- data$EntrezID[data[[proxyScore]] == "HighConf"]
  
  nonHits <- setdiff(data$EntrezID, Hits)
  
  outPrefix <- paste("KEGG", iteration, sep = "_")
  
  # 1) Contraction - [Pathway Analysis]
  message(getwd())
  data <- ComputeEnrichment(pathwayData, Hits, nonHits, outPrefix, data, iteration)
  data <- data.frame(data, temp = 0, stringsAsFactors = FALSE)
  kName1 <- paste0("KEGG.class.iteration", iteration)
  kName2 <- paste0("KEGG.", iteration)
  names(data)[names(data) == "temp"] <- kName1
  data[[kName1]][data$KEGG == "Yes" & (data[[proxyScore]] %in% c("MedConf", "HighConf"))] <- "HighConf"
  data[[kName1]][data$KEGG != "Yes" & (data[[proxyScore]] %in% c("MedConf", "HighConf"))] <- "MedConf"
  names(data)[names(data) == "KEGG"] <- kName2
  
  hit.Genes <- data$EntrezID[data[[kName1]] == "HighConf"]
  myOrignalGenes <- data$GeneSymbol[data[[kName1]] == "HighConf"]
  
  # 2) Expansion - [Network Analysis]
  message("*", paste0(triage.wd, "Network_iteration_V3_longform.R"), "**")
  
  source(paste0(triage.wd, "Network_iteration_V3.R"), local = TRUE)
  data <- data.frame(data, temp1 = "No", temp2 = data[[kName1]], stringsAsFactors = FALSE)
  nName1 <- paste0("Network.", iteration)
  nName2 <- paste0("Network.class.iteration", iteration)
  names(data)[names(data) == "temp1"] <- nName1
  names(data)[names(data) == "temp2"] <- nName2
  data[[nName1]][data$EntrezID %in% gNames2] <- "Yes"
  # data[[nName2]][data$EntrezID %in% gNames2 & data[[kName1]] > 0] <- 1
  data[[nName2]][data$EntrezID %in% gNames2 & data[[kName1]] %in% c("MedConf", "HighConf")] <- "HighConf"
  
  
  
  
  
  if((iteration != 1 && identical(data[[nName2]], data[[paste0("Network.class.iteration", iteration-1)]])) 
     || (iteration >= 5 
         && identical(data[[nName2]], data[[paste0("Network.class.iteration", iteration-2)]]) 
         && identical(data[[paste0("Network.class.iteration", iteration-1)]], data[[paste0("Network.class.iteration", iteration-3)]])
         && (length(data$EntrezID[data[[nName2]]== "HighConf"]) > length(data$EntrezID[data[[paste0("Network.class.iteration", iteration-1)]] == "HighConf"])))) { 
    #dupCols <- (ncol(data)-1):ncol(data)      #This just shaves off the last two columns 
    #data <- data[, -dupCols]
    counter <- FALSE
  }
  
  #print(paste("iteration: ", iteration, "\n"))
  message(paste("iteration: ", iteration, "\n"))
  proxyScore <- nName2
  iteration <- iteration + 1
} # end of while loop

### Append Enrichment Info --------------------------------------------
# set the cutoff for the enriched pathways to display
pval_threshold <- 0.05

iterationNum <<- iteration - 1
enrichFileName <- paste0(outPrefix,".Enrichment_", iterationNum, ".csv")
pathEnrich <- read.csv(enrichFileName, stringsAsFactors = FALSE)


#pathEnrich <- pathEnrich[pathEnrich$pVal < pval_threshold, ]  #Removed so that all pathways are listed even those that don't pass significance test

tempL <- strsplit(pathEnrich$HitGeneNames, split = ", ")
names(tempL) <- pathEnrich$Pathway

library(reshape2)
tempL <- lapply(seq(tempL), function(i) {
  m <- melt(tempL[i])
  names(m) <- c("GeneSymbol", paste0("Pathway", i))
  return(m)
})

tempDF <- Reduce(function(x, y) merge(x, y, by = "GeneSymbol", all = T), tempL)
#samp <<- tempDF[, 2:ncol(tempDF)]
samp <- as.data.frame(tempDF[, -1])
pathVector <- sapply(seq(nrow(samp)), function(i) unlist(paste(samp[i, which(!is.na(samp[i, ]))], collapse = ", "))) #Switched to comma seperator from " ; " 

pathDF <- data.frame(GeneSymbol = tempDF$GeneSymbol, Pathway = pathVector, stringsAsFactors = F)

## Save the results into output files in the TRIAGEoutputFiles folder
out <- merge(data, pathDF, by = "GeneSymbol", all = T)


message(getwd(), "#####")

######## Define names and write final files
outputFileName <- paste0(AnalysisTitle, "_TRIAGEoutput_ALL.csv")
write.csv(out, file = outputFileName, row.names = F)
final_enriched_pathway_file <- paste0(AnalysisTitle, "_TRIAGE_enrichment_final", ".csv")
write.csv(pathEnrich, file = final_enriched_pathway_file, row.names = F)

#Create TRIAGE output
TRIAGEoutput <- read.csv(outputFileName, stringsAsFactors = F)

FinalIterationNetworkColumn <- paste0("Network.class.iteration", iterationNum)

TRIAGEoutput <- TRIAGEoutput %>%
  mutate(TRIAGEhit = ifelse(get(FinalIterationNetworkColumn, envir = as.environment(TRIAGEoutput)) == "HighConf", 
                            "Yes",
                            ""))

########################
######## Pathway Output
#########################
FinalEnrichment.df <- pathEnrich

TRIAGEhits.highConf.matrix.GS <- matrix(TRIAGEoutput$GeneSymbol[TRIAGEoutput$TRIAGEhit == "Yes" & TRIAGEoutput$ConfidenceCategory == "HighConf"])
TRIAGEhits.medConf.matrix.GS <- matrix(TRIAGEoutput$GeneSymbol[TRIAGEoutput$TRIAGEhit == "Yes" & TRIAGEoutput$ConfidenceCategory == "MedConf"])


#Genrate columns with high confidence and med confidence (based on input) genes of each pathway.
FinalEnrichment.df$HighScoreGenes <- NA
FinalEnrichment.df$HighScoreGenesNames <- NA
FinalEnrichment.df$MedScoreGenesNames <- NA

for (i in 1:length(FinalEnrichment.df$Genes)) {
  temp.path.string <- unlist(strsplit(FinalEnrichment.df$HitGeneNames[i], ", "))
  out.HC <- paste(intersect(temp.path.string, TRIAGEhits.highConf.matrix.GS),collapse = ", ")
  out.HC.count <- length(intersect(temp.path.string, TRIAGEhits.highConf.matrix.GS))
  out.MC <- paste(intersect(temp.path.string, TRIAGEhits.medConf.matrix.GS),collapse = ", ")
  FinalEnrichment.df$HighScoreGenes[i] <- out.HC.count
  FinalEnrichment.df$HighScoreGenesNames[i] <- out.HC
  FinalEnrichment.df$MedScoreGenesNames[i] <- out.MC
}


########### Generate Enrichment Score

FinalEnrichment.df$EnrichScore <- NA

for (i in 1:length(FinalEnrichment.df$Pathway)) {
  GeneHitGeneRatio <- FinalEnrichment.df$HitGenes[i] / FinalEnrichment.df$Genes[i]
  HighConfHitGeneRation <-  FinalEnrichment.df$HighScoreGenes[i] / FinalEnrichment.df$HitGenes[i]
  FinalEnrichment.df$EnrichScore[i] <- round(((GeneHitGeneRatio + HighConfHitGeneRation) / 2), 3)
}


FinalEnrichment.condensed <- FinalEnrichment.df[, c("Pathway", "pVal", "pValFDR", "pValBonferroni", "Genes", "HitGenes", "HighScoreGenes", "HighScoreGenesNames", "MedScoreGenesNames", "EnrichScore")]
############# Write files to new Directory
TRIAGEoutput.dir <- paste0(parnas.wd, "TRIAGEoutput/" )
setwd(TRIAGEoutput.dir)

TRIAGE.cond.output.name <- paste0(AnalysisTitle, "_", "TRIAGE_hits.csv")
Enrichment.cond.output.name <- paste0(AnalysisTitle, "_", "TRIAGE_enrichment.csv")

write.csv(TRIAGEoutput, file = TRIAGE.cond.output.name)
write.csv(FinalEnrichment.condensed, file = Enrichment.cond.output.name)
  
######################--
### Compare pathways of top score versus triage
######################--
setwd(paste0(parnas.wd, "TRIAGEoutput/"))
  
TopScore_path.df <- read.csv("parnas_crispr_ZscorePathwaysAnalysis.csv", stringsAsFactors = F)
TRIAGE_path.df <- read.csv("parnas_crispr_TRIAGE_enrichment.csv", stringsAsFactors = F)


#Get matrices
TopScore_path.matrix <- TopScore_path.df$Pathway[which(TopScore_path.df$pValFDR <= 0.1)]
TRIAGE_path.matrix <- TRIAGE_path.df$Pathway[which(TRIAGE_path.df$pValFDR <= 0.1)]

#Get TRIAGE unique pathways
TRIAGE_path_unique <- setdiff(TRIAGE_path.matrix, TopScore_path.matrix)


################---
##### Get TRIAGE Hits ----
###############---


setwd(TRIAGEoutput.dir)





## Clear TRIAGE data frame for critical columns only

condensed.columns <- c("GeneSymbol", "EntrezID", "log2FoldChange", "pvalue", "Rank", "TRIAGEscore", "ConfidenceCategory", "Pathway", "TRIAGEhit")


parnas.triage <- read.csv("parnas_crispr_TRIAGE_hits.csv", stringsAsFactors = F)
  
parnas.triage <- parnas.triage[, condensed.columns]
  

##################---
# True Positive and true negative ----
#################---

validation.lower.bound <- 0.5
validation.upper.bound <- 2.5



# Distribution of scores for secondary screen
parnas.sec.histo <- ggplot(parnas.secondary.df)+
  geom_histogram(aes(x = zscore), binwidth = .5, 
                 col="black", 
                 size=.5, fill="white")+
  ylab("Frequency") + xlab("Zscore")+
  geom_rect(aes(xmin=0.5, xmax=3, ymin=-10, ymax=1100), fill = "white", color="red", alpha=0, size = 0.8)+
  #geom_vline(xintercept = c(0.7, 3), color = "red")+
  ggtitle(paste0("Parnase et al, Secondary Screen", ": Score Distribution"))+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=10),
        axis.title=element_text(size=20,face="bold"),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.line.x = element_line(colour = 'black', size = 2),
        panel.background = element_blank())

setwd(paste0(parnas.wd, "Figures/"))
ggsave("Histogram_secondary_parnas.png", width = 12.5, height = 12.5)


####Secondary Screen - Pub figure ----
parnas.sec.histo <- ggplot(parnas.secondary.df)+
  geom_histogram(aes(x = zscore), binwidth = .25, 
                 col="white", 
                 size=.5, fill="black")+
  ylab("Frequency") + xlab("Log 2 Fold Change")+
  #geom_rect(aes(xmin=validation.lower.bound, xmax=validation.upper.bound, ymin=-10, ymax=700), fill = "white", color="red", alpha=0, size = 2, linejoin = "round", linetype = "dashed")+
  geom_vline(xintercept = c(validation.lower.bound, validation.upper.bound), color = "red", size = 2, linetype = "dashed")+
  ggtitle(paste0("Parnase et al, Secondary Screen", ": Score Distribution"))+
  scale_x_continuous(breaks = c(-2, -1, 0, validation.lower.bound, 1, 2, validation.upper.bound),
                     labels = c(-2, -1, 0, validation.lower.bound, 1, 2, validation.upper.bound), limits = c(-3, 3))+
  theme_tufte(base_size = 20)+
  theme(legend.position="none", axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20, colour = "black", angle = 0, hjust = 0.5),
        axis.title=element_text(size=40,face="bold"),
        axis.text.y.left = element_text(size=30),
        axis.title.x = element_text(vjust = -1),
        axis.line.y = element_line(colour = 'black', size = 2, lineend = "round"),
        axis.line.x = element_line(colour = 'black', size = 2, lineend = "round"),
        strip.text.x = element_text(size = 25),
        aspect.ratio = 1/1,
        panel.background = element_blank())


print(parnas.sec.histo)

setwd(paste0(parnas.wd, "Figures/"))
ggsave("Secondary_screen_parnas_manuscript.png", width = 12.5, height = 12.5)



###############--
# Generate calssification tabs ----
##############--

validation.test.count <- 100 # NUmber of cutoffs to be chosen for validation test



validation.cutoffs <- seq(validation.lower.bound, validation.upper.bound, length.out = validation.test.count)

analysis.types <- c("all.hits"
                    , "triage.hits"
                    )



analysis.labels <- c("Selected Hits", "TRIAGE Hits")


#CHaractarize membership in original input

parnas.df.membership <- parnas.triage %>%
  mutate(top.hits = ifelse(TRIAGEscore == 1, "hit", "nonhit")) %>%
  mutate(all.hits = ifelse(TRIAGEscore >= 0.5 , "hit", "nonhit")) %>%
  mutate(triage.hits = ifelse(TRIAGEhit == "Yes", "hit", "nonhit")) %>%
  dplyr::select(GeneSymbol, EntrezID, log2FoldChange, pvalue, top.hits, all.hits, triage.hits)





### Create dataframe of values

prediction.columns <- c("hit.selection", "cutoff", "tp", "fp", "tn", "fn", "prev", "pred_pos", "tpr", "fpr", "acc", "err", "fnr", "spec", "prec", "npv", "plr", "nlr", "dor", "youdens", "f1", "fdr", "f_or")
prediction.colnames <- c("Hit Selection", "Validation Cutoff", "True Positives", "False Positives", "True Negatives", "False Negatives", "Prevelance", "Predicted Positive", "True Positive Rate"
                       , "False Positive Rate", "Accuracy", "Error", "False Negative Rate", "Specificity", "Precision", "Negative Predictive Value", "Positive Likelihood Ratio", "Negative Likelihood Ratio"
                       , "Diagnostic Odds Ration (DOR)", "Youden's J statistic", "F1 Score ", "False Discovery Rate (FDR)", "False Ommision Rate (FOR)")

prediction.df <- data.frame(matrix(ncol = length(prediction.columns), nrow = 0))
colnames(prediction.df) <- prediction.columns

for (a in 1:length(analysis.types)) {
  analysis.type <- analysis.types[a]
  
  for (c in 1:validation.test.count){      #generate a matrix tab o
    validation.cutoff <- validation.cutoffs[c]
    
    ### Add a row in the prediction data frame
    prediction.df[nrow(prediction.df)+1, ] <- NA
    
    #### fill in analysis type and cutoff
    prediction.df$hit.selection[nrow(prediction.df)] <- analysis.type
    
    prediction.df$cutoff[nrow(prediction.df)] <- validation.cutoff
    
    #Find genes that are validated
    secondary.screen <- parnas.secondary.df %>%
      mutate(cutoff.val = ifelse(zscore >= validation.cutoff, "hit", "nonhit"))
    
    validated.genesymbol.matrix <- secondary.screen$GeneSymbol[which(secondary.screen$cutoff.val == "hit")]
    
    #Assign validation call to input file
    primary.screen <- parnas.df.membership %>%
      mutate(validated = ifelse(GeneSymbol %in% validated.genesymbol.matrix, "hit", "nonhit"))
    
    #Create dataframe
    
    lvs <- c("nonhit", "hit")    
    truth <- factor(primary.screen$validated,
                    levels = rev(lvs))
    pred <- factor(primary.screen[[analysis.type]],
      levels = rev(lvs))
    
    xtab <- table(pred, truth)
    
    #Create confusion matrix
    
    error.matrix <- confusionMatrix(xtab)
    
    
    #
    
    
    ###### Fill in all the prediction values
    
    #True positive
    prediction.df$tp[nrow(prediction.df)] <- error.matrix$table["hit", "hit"]
    
    #True negative
    prediction.df$tn[nrow(prediction.df)] <- error.matrix$table["nonhit", "nonhit"]
    
    #False positive
    prediction.df$fp[nrow(prediction.df)] <- error.matrix$table["hit", "nonhit"]
    
    #False negative
    prediction.df$fn[nrow(prediction.df)] <- error.matrix$table["nonhit", "hit"]
    
    #Prevelance
    prediction.df$prev[nrow(prediction.df)] <- (error.matrix$table["hit", "hit"] + error.matrix$table["nonhit", "hit"]) / sum(error.matrix$table)
    
    #Predicted Positive
    prediction.df$pred_pos[nrow(prediction.df)] <- error.matrix$table["hit", "hit"] + error.matrix$table["hit", "nonhit"]
    
    
    # True Positive rate
    prediction.df$tpr[nrow(prediction.df)] <- error.matrix$byClass["Sensitivity"]
    
    #False positive rate 
    prediction.df$fpr[nrow(prediction.df)] <- 1 - error.matrix$byClass["Specificity"]
    
    # Accuracy 
    prediction.df$acc[nrow(prediction.df)] <- error.matrix$overall["Accuracy"]
    
    # Error
    prediction.df$err[nrow(prediction.df)] <- 1 - error.matrix$overall["Accuracy"]
    
    
    # False negative rate
    prediction.df$fnr[nrow(prediction.df)] <- 1 - error.matrix$byClass["Sensitivity"]
    
    #Specificity
    prediction.df$spec[nrow(prediction.df)] <- error.matrix$byClass["Specificity"]
    
    #Precision (PPV)
    prediction.df$prec[nrow(prediction.df)] <- error.matrix$byClass["Precision"]
    
    # Negative predictive value 
    prediction.df$npv[nrow(prediction.df)] <- error.matrix$byClass["Neg Pred Value"]
    
    # Positive Likelihood ration
    prediction.df$plr[nrow(prediction.df)] <- error.matrix$byClass["Sensitivity"] / (1 - error.matrix$byClass["Specificity"])
    
    # Negative likelihood ratio
    prediction.df$nlr[nrow(prediction.df)] <- (1 - error.matrix$byClass["Sensitivity"]) / error.matrix$byClass["Specificity"]
    
    # Diagnostic odds ratio
    prediction.df$dor[nrow(prediction.df)] <- (error.matrix$byClass["Sensitivity"] / (1 - error.matrix$byClass["Specificity"])) / ((1 - error.matrix$byClass["Sensitivity"]) / error.matrix$byClass["Specificity"])
    
    # Youden's J statistic
    prediction.df$youdens[nrow(prediction.df)] <- error.matrix$byClass["Sensitivity"] + error.matrix$byClass["Specificity"] - 1
    
    # F1 score
    prediction.df$f1[nrow(prediction.df)] <- 2 * ((error.matrix$byClass["Precision"] * error.matrix$byClass["Sensitivity"]) / (error.matrix$byClass["Precision"] + error.matrix$byClass["Sensitivity"]))
    
    # FDR
    prediction.df$fdr[nrow(prediction.df)] <- (error.matrix$table["hit", "nonhit"])/ (error.matrix$table["hit", "hit"] + error.matrix$table["hit", "nonhit"])
    
    #FOR
    prediction.df$f_or[nrow(prediction.df)] <- (error.matrix$table["nonhit", "hit"])/ (error.matrix$table["nonhit", "hit"] + error.matrix$table["nonhit", "nonhit"])
    
  }
}



############--
## Prediction figures, A ----
###########--

prediction.measures <- prediction.columns[!prediction.columns %in% c("hit.selection")]
prediction.labels <- prediction.colnames[!prediction.colnames %in% c("Hit Selection")]




for (x in 1:(length(prediction.measures)/2)){
  prediction.measure.x <- prediction.measures[x]
  prediction.label.x <- prediction.labels[x]
  
  prediction.measures.y <- prediction.measures[!prediction.measures %in% prediction.measure.x]
  prediction.labels.y <- prediction.labels[!prediction.labels %in% prediction.label.x]
  for (y in 1:length(prediction.measures.y)){
    prediction.measure.y <- prediction.measures.y[y]
    prediction.label.y <- prediction.labels[y]
    
    ## FIgure name
    figure.name <- paste0(prediction.measure.y, "_vs_", prediction.measure.x)
    figure.title <- paste0(prediction.label.y, " vs. ", prediction.label.x)
    
    
    #Set aesthetcis
    hit.selection.colors <- c("#D5AC4C", "black",  "#1A4384")
    hit.selection.colors2 <- c("midnightblue", "black",  "firebrick")
    
    
    #Generate figure
    
    prediction.figure <- ggplot(prediction.df, aes(get(prediction.measure.x), get(prediction.measure.y), color=factor(hit.selection)))+
      geom_line(size = 1, alpha = 0.6) +
      geom_point(size = 2)+
      # scale_y_continuous(breaks=c(min(BROAD.tlr.plus.entrez.ordered$ZHL),cutoff,0,max(BROAD.tlr.plus.entrez.ordered$ZHL)),
      #                    labels =round(c(min(BROAD.tlr.plus.entrez.ordered$ZHL),cutoff,0,max(BROAD.tlr.plus.entrez.ordered$ZHL)),1) ) +
      # scale_x_continuous(breaks = c(1, nrow(BROAD.tlr.plus.entrez.ordered) ))+
      scale_color_manual(values = hit.selection.colors2, labels = analysis.labels)+
      labs(title = paste0("CRISPR/Cas9 TNF Validation Screen    -     ", figure.title), y = prediction.label.y, x = prediction.label.x, color = "Hit Selection Method")+
      # annotate(geom="text", x=0.1*nrow(BROAD.tlr.plus.entrez.ordered), y=0.75*cutoff, label="2.5% cutoff",
      #          color="black", size =7)+
      # annotate(geom="text", x=0.1*nrow(BROAD.tlr.plus.entrez.ordered), y=1.25*cutoff, 
      #          label=paste0(length(which(BROAD.tlr.plus.entrez.ordered$ZHL >= cutoff)), " hits"),
      #          color="black", size =7)+
      theme(axis.text.y = element_text(size=15),
            axis.text.x = element_text(size=15),
            axis.title=element_text(size=20,face="bold"),
            axis.line.y = element_line(colour = 'black', size = 2),
            axis.line.x = element_line(colour = 'black', size = 2),
            panel.background = element_blank(),   
            panel.grid.major = element_line(colour = "grey"), 
            panel.grid.minor = element_line(colour = "grey", size = 0.25),
            plot.title = element_text(size = 20, hjust = 0.5, vjust = 2, 
                                      margin = margin(b = 6, t = 6)))
    
    setwd(paste0(parnas.wd, "Figures/PredictionFigures/AbsolutePrediction_CompleteScreen/"))
    ggsave(paste0(figure.name, "_AbsPred_CompScreen.png"), width = 13.5, height = 12.5)
  }
  
}



############### --
# Prediction Figures -Delta ----
############## --- 


############# ---
#Prediction Figures - versus validation range----
#############----
validation.test.count <- 100 # NUmber of cutoffs to be chosen for validation test

validation.lower.bound <- 0.5
validation.upper.bound <- 2.5

validation.cutoffs <- seq(validation.lower.bound, validation.upper.bound, length.out = validation.test.count)

analysis.types <- c("all.hits"
                    , "triage.hits"
)



analysis.labels <- c("Selected Hits", "TRIAGE Hits")


#CHaractarize membership in original input

parnas.df.membership <- parnas.triage %>%
  mutate(top.hits = ifelse(TRIAGEscore == 1, "hit", "nonhit")) %>%
  mutate(all.hits = ifelse(TRIAGEscore >= 0.5 , "hit", "nonhit")) %>%
  mutate(triage.hits = ifelse(TRIAGEhit == "Yes", "hit", "nonhit")) %>%
  dplyr::select(GeneSymbol, EntrezID, log2FoldChange, pvalue, top.hits, all.hits, triage.hits)





### Create dataframe of values

prediction.columns <- c("hit.selection", "cutoff", "tp", "fp", "tn", "fn", "prev", "pred_pos", "tpr", "fpr", "acc", "err", "fnr", "spec", "prec", "npv", "plr", "nlr", "dor", "youdens", "f1", "fdr", "f_or")
prediction.colnames <- c("Hit Selection", "Validation Cutoff", "True Positives", "False Positives", "True Negatives", "False Negatives", "Prevelance", "Predicted Positive", "True Positive Rate"
                         , "False Positive Rate", "Accuracy", "Error", "False Negative Rate", "Specificity", "Precision", "Negative Predictive Value", "Positive Likelihood Ratio", "Negative Likelihood Ratio"
                         , "Diagnostic Odds Ration (DOR)", "Youden's J statistic", "F1 Score ", "False Discovery Rate (FDR)", "False Ommision Rate (FOR)")

prediction.df <- data.frame(matrix(ncol = length(prediction.columns), nrow = 0))
colnames(prediction.df) <- prediction.columns

for (a in 1:length(analysis.types)) {
  analysis.type <- analysis.types[a]
  
  for (c in 1:validation.test.count){      #generate a matrix tab o
    validation.cutoff <- validation.cutoffs[c]
    
    ### Add a row in the prediction data frame
    prediction.df[nrow(prediction.df)+1, ] <- NA
    
    #### fill in analysis type and cutoff
    prediction.df$hit.selection[nrow(prediction.df)] <- analysis.type
    
    prediction.df$cutoff[nrow(prediction.df)] <- validation.cutoff
    
    #Find genes that are validated
    secondary.screen <- parnas.secondary.df %>%
      mutate(cutoff.val = ifelse(zscore >= validation.cutoff, "hit", "nonhit"))
    
    validated.genesymbol.matrix <- secondary.screen$GeneSymbol[which(secondary.screen$cutoff.val == "hit")]
    
    #Assign validation call to input file
    primary.screen <- parnas.df.membership %>%
      mutate(validated = ifelse(GeneSymbol %in% validated.genesymbol.matrix, "hit", "nonhit"))
    
    #Create dataframe
    
    lvs <- c("nonhit", "hit")    
    truth <- factor(primary.screen$validated,
                    levels = rev(lvs))
    pred <- factor(primary.screen[[analysis.type]],
                   levels = rev(lvs))
    
    xtab <- table(pred, truth)
    
    #Create confusion matrix
    
    error.matrix <- confusionMatrix(xtab)
    
    
    #
    
    
    ###### Fill in all the prediction values
    
    #True positive
    prediction.df$tp[nrow(prediction.df)] <- error.matrix$table["hit", "hit"]
    
    #True negative
    prediction.df$tn[nrow(prediction.df)] <- error.matrix$table["nonhit", "nonhit"]
    
    #False positive
    prediction.df$fp[nrow(prediction.df)] <- error.matrix$table["hit", "nonhit"]
    
    #False negative
    prediction.df$fn[nrow(prediction.df)] <- error.matrix$table["nonhit", "hit"]
    
    #Prevelance
    prediction.df$prev[nrow(prediction.df)] <- (error.matrix$table["hit", "hit"] + error.matrix$table["nonhit", "hit"]) / sum(error.matrix$table)
    
    #Predicted Positive
    prediction.df$pred_pos[nrow(prediction.df)] <- error.matrix$table["hit", "hit"] + error.matrix$table["hit", "nonhit"]
    
    
    # True Positive rate
    prediction.df$tpr[nrow(prediction.df)] <- error.matrix$byClass["Sensitivity"]
    
    #False positive rate 
    prediction.df$fpr[nrow(prediction.df)] <- 1 - error.matrix$byClass["Specificity"]
    
    # Accuracy 
    prediction.df$acc[nrow(prediction.df)] <- error.matrix$overall["Accuracy"]
    
    # Error
    prediction.df$err[nrow(prediction.df)] <- 1 - error.matrix$overall["Accuracy"]
    
    
    # False negative rate
    prediction.df$fnr[nrow(prediction.df)] <- 1 - error.matrix$byClass["Sensitivity"]
    
    #Specificity
    prediction.df$spec[nrow(prediction.df)] <- error.matrix$byClass["Specificity"]
    
    #Precision (PPV)
    prediction.df$prec[nrow(prediction.df)] <- error.matrix$byClass["Precision"]
    
    # Negative predictive value 
    prediction.df$npv[nrow(prediction.df)] <- error.matrix$byClass["Neg Pred Value"]
    
    # Positive Likelihood ration
    prediction.df$plr[nrow(prediction.df)] <- error.matrix$byClass["Sensitivity"] / (1 - error.matrix$byClass["Specificity"])
    
    # Negative likelihood ratio
    prediction.df$nlr[nrow(prediction.df)] <- (1 - error.matrix$byClass["Sensitivity"]) / error.matrix$byClass["Specificity"]
    
    # Diagnostic odds ratio
    prediction.df$dor[nrow(prediction.df)] <- (error.matrix$byClass["Sensitivity"] / (1 - error.matrix$byClass["Specificity"])) / ((1 - error.matrix$byClass["Sensitivity"]) / error.matrix$byClass["Specificity"])
    
    # Youden's J statistic
    prediction.df$youdens[nrow(prediction.df)] <- error.matrix$byClass["Sensitivity"] + error.matrix$byClass["Specificity"] - 1
    
    # F1 score
    prediction.df$f1[nrow(prediction.df)] <- 2 * ((error.matrix$byClass["Precision"] * error.matrix$byClass["Sensitivity"]) / (error.matrix$byClass["Precision"] + error.matrix$byClass["Sensitivity"]))
    
    # FDR
    prediction.df$fdr[nrow(prediction.df)] <- (error.matrix$table["hit", "nonhit"])/ (error.matrix$table["hit", "hit"] + error.matrix$table["hit", "nonhit"])
    
    #FOR
    prediction.df$f_or[nrow(prediction.df)] <- (error.matrix$table["nonhit", "hit"])/ (error.matrix$table["nonhit", "hit"] + error.matrix$table["nonhit", "nonhit"])
    
  }
}



############--
## Prediction figures, A 
###########--

prediction.measures <- prediction.columns[!prediction.columns %in% c("hit.selection", "cutoff")]
prediction.labels <- prediction.colnames[!prediction.colnames %in% c("Hit Selection", "Validation Cutoff")]




for (x in 1:length(prediction.measures)){
  prediction.measure.x <- prediction.measures[x]
  prediction.label.x <- prediction.labels[x]
  

    
    ## FIgure name
    figure.name <- paste0(prediction.measure.x, "_vs_cutoff")
    figure.title <- paste0(prediction.label.x, " vs. Cuotff")
    
    
    #Set aesthetcis
    hit.selection.colors <- c("black",  "#1A4384")
    hit.selection.colors2 <- c("black",  "firebrick")
    
    
    #Generate figure
    
    prediction.figure <- ggplot(prediction.df, aes(get(prediction.measure.x), cutoff, color=factor(hit.selection)))+
      geom_line(size = 1, alpha = 0.6) +
      geom_point(size = 2)+
      # scale_y_continuous(breaks=c(min(BROAD.tlr.plus.entrez.ordered$ZHL),cutoff,0,max(BROAD.tlr.plus.entrez.ordered$ZHL)),
      #                    labels =round(c(min(BROAD.tlr.plus.entrez.ordered$ZHL),cutoff,0,max(BROAD.tlr.plus.entrez.ordered$ZHL)),1) ) +
      # scale_x_continuous(breaks = c(1, nrow(BROAD.tlr.plus.entrez.ordered) ))+
      scale_color_manual(values = hit.selection.colors2, labels = analysis.labels)+
      labs(title = paste0("CRISPR/Cas9 TNF Validation Screen    -     ", figure.title), y = "Validation Cutoffs", x = prediction.label.x, color = "Hit Selection Method")+
      # annotate(geom="text", x=0.1*nrow(BROAD.tlr.plus.entrez.ordered), y=0.75*cutoff, label="2.5% cutoff",
      #          color="black", size =7)+
      # annotate(geom="text", x=0.1*nrow(BROAD.tlr.plus.entrez.ordered), y=1.25*cutoff, 
      #          label=paste0(length(which(BROAD.tlr.plus.entrez.ordered$ZHL >= cutoff)), " hits"),
      #          color="black", size =7)+
      theme(axis.text.y = element_text(size=15),
            axis.text.x = element_text(size=15),
            axis.title=element_text(size=20,face="bold"),
            axis.line.y = element_line(colour = 'black', size = 2),
            axis.line.x = element_line(colour = 'black', size = 2),
            panel.background = element_blank(),   
            panel.grid.major = element_line(colour = "grey"), 
            panel.grid.minor = element_line(colour = "grey", size = 0.25),
            plot.title = element_text(size = 20, hjust = 0.5, vjust = 2, 
                                      margin = margin(b = 6, t = 6)))
    
    setwd(paste0(parnas.wd, "Figures/PredictionFigures/AbsolutePrediction_versusCutoff/"))
    ggsave(paste0(figure.name, "_ParnasSecondary.png"), width = 13.5, height = 12.5)
}

  

############# ---
#Prediction Figures - versus validation range, input secondary only----
#############----
validation.test.count <- 20 # NUmber of cutoffs to be chosen for validation test

validation.lower.bound <- 0.5
validation.upper.bound <- 2.5

validation.cutoffs <- seq(validation.lower.bound, validation.upper.bound, length.out = validation.test.count)

analysis.types <- c("all.hits"
                    , "triage.hits"
)



analysis.labels <- c("Selected Hits", "TRIAGE Hits")


#CHaractarize membership in original input

parnas.df.membership <- parnas.triage %>%
  mutate(top.hits = ifelse(TRIAGEscore == 1, "hit", "nonhit")) %>%
  mutate(all.hits = ifelse(TRIAGEscore >= 0.5 , "hit", "nonhit")) %>%
  mutate(triage.hits = ifelse(TRIAGEhit == "Yes", "hit", "nonhit")) %>%
  filter(GeneSymbol %in% parnas.secondary.df$GeneSymbol) %>%
  dplyr::select(GeneSymbol, EntrezID, log2FoldChange, pvalue, top.hits, all.hits, triage.hits)





### Create dataframe of values

prediction.columns <- c("hit.selection", "cutoff", "tp", "fp", "tn", "fn", "prev", "pred_pos", "tpr", "fpr", "acc", "err", "fnr", "spec", "prec", "npv", "plr", "nlr", "dor", "youdens", "f1", "fdr", "f_or")
prediction.colnames <- c("Hit Selection", "Validation Cutoff", "True Positives", "False Positives", "True Negatives", "False Negatives", "Prevelance", "Predicted Positive", "True Positive Rate"
                         , "False Positive Rate", "Accuracy", "Error", "False Negative Rate", "Specificity", "Precision", "Negative Predictive Value", "Positive Likelihood Ratio", "Negative Likelihood Ratio"
                         , "Diagnostic Odds Ration (DOR)", "Youden's J statistic", "F1 Score ", "False Discovery Rate (FDR)", "False Ommision Rate (FOR)")

prediction.df <- data.frame(matrix(ncol = length(prediction.columns), nrow = 0))
colnames(prediction.df) <- prediction.columns

for (a in 1:length(analysis.types)) {
  analysis.type <- analysis.types[a]
  
  for (c in 1:validation.test.count){      #generate a matrix tab o
    validation.cutoff <- validation.cutoffs[c]
    
    ### Add a row in the prediction data frame
    prediction.df[nrow(prediction.df)+1, ] <- NA
    
    #### fill in analysis type and cutoff
    prediction.df$hit.selection[nrow(prediction.df)] <- analysis.type
    
    prediction.df$cutoff[nrow(prediction.df)] <- validation.cutoff
    
    #Find genes that are validated
    secondary.screen <- parnas.secondary.df %>%
      mutate(cutoff.val = ifelse(zscore >= validation.cutoff, "hit", "nonhit"))
    
    validated.genesymbol.matrix <- secondary.screen$GeneSymbol[which(secondary.screen$cutoff.val == "hit")]
    
    #Assign validation call to input file
    primary.screen <- parnas.df.membership %>%
      mutate(validated = ifelse(GeneSymbol %in% validated.genesymbol.matrix, "hit", "nonhit"))
    
    #Create dataframe
    
    lvs <- c("nonhit", "hit")    
    truth <- factor(primary.screen$validated,
                    levels = rev(lvs))
    pred <- factor(primary.screen[[analysis.type]],
                   levels = rev(lvs))
    
    xtab <- table(pred, truth)
    
    #Create confusion matrix
    
    error.matrix <- confusionMatrix(xtab)
    
    
    #
    
    
    ###### Fill in all the prediction values
    
    #True positive
    prediction.df$tp[nrow(prediction.df)] <- error.matrix$table["hit", "hit"]
    
    #True negative
    prediction.df$tn[nrow(prediction.df)] <- error.matrix$table["nonhit", "nonhit"]
    
    #False positive
    prediction.df$fp[nrow(prediction.df)] <- error.matrix$table["hit", "nonhit"]
    
    #False negative
    prediction.df$fn[nrow(prediction.df)] <- error.matrix$table["nonhit", "hit"]
    
    #Prevelance
    prediction.df$prev[nrow(prediction.df)] <- (error.matrix$table["hit", "hit"] + error.matrix$table["nonhit", "hit"]) / sum(error.matrix$table)
    
    #Predicted Positive
    prediction.df$pred_pos[nrow(prediction.df)] <- error.matrix$table["hit", "hit"] + error.matrix$table["hit", "nonhit"]
    
    
    # True Positive rate
    prediction.df$tpr[nrow(prediction.df)] <- error.matrix$byClass["Sensitivity"]
    
    #False positive rate 
    prediction.df$fpr[nrow(prediction.df)] <- 1 - error.matrix$byClass["Specificity"]
    
    # Accuracy 
    prediction.df$acc[nrow(prediction.df)] <- error.matrix$overall["Accuracy"]
    
    # Error
    prediction.df$err[nrow(prediction.df)] <- 1 - error.matrix$overall["Accuracy"]
    
    
    # False negative rate
    prediction.df$fnr[nrow(prediction.df)] <- 1 - error.matrix$byClass["Sensitivity"]
    
    #Specificity
    prediction.df$spec[nrow(prediction.df)] <- error.matrix$byClass["Specificity"]
    
    #Precision (PPV)
    prediction.df$prec[nrow(prediction.df)] <- error.matrix$byClass["Precision"]
    
    # Negative predictive value 
    prediction.df$npv[nrow(prediction.df)] <- error.matrix$byClass["Neg Pred Value"]
    
    # Positive Likelihood ration
    prediction.df$plr[nrow(prediction.df)] <- error.matrix$byClass["Sensitivity"] / (1 - error.matrix$byClass["Specificity"])
    
    # Negative likelihood ratio
    prediction.df$nlr[nrow(prediction.df)] <- (1 - error.matrix$byClass["Sensitivity"]) / error.matrix$byClass["Specificity"]
    
    # Diagnostic odds ratio
    prediction.df$dor[nrow(prediction.df)] <- (error.matrix$byClass["Sensitivity"] / (1 - error.matrix$byClass["Specificity"])) / ((1 - error.matrix$byClass["Sensitivity"]) / error.matrix$byClass["Specificity"])
    
    # Youden's J statistic
    prediction.df$youdens[nrow(prediction.df)] <- error.matrix$byClass["Sensitivity"] + error.matrix$byClass["Specificity"] - 1
    
    # F1 score
    prediction.df$f1[nrow(prediction.df)] <- 2 * ((error.matrix$byClass["Precision"] * error.matrix$byClass["Sensitivity"]) / (error.matrix$byClass["Precision"] + error.matrix$byClass["Sensitivity"]))
    
    # FDR
    prediction.df$fdr[nrow(prediction.df)] <- (error.matrix$table["hit", "nonhit"])/ (error.matrix$table["hit", "hit"] + error.matrix$table["hit", "nonhit"])
    
    #FOR
    prediction.df$f_or[nrow(prediction.df)] <- (error.matrix$table["nonhit", "hit"])/ (error.matrix$table["nonhit", "hit"] + error.matrix$table["nonhit", "nonhit"])
    
  }
}



############--
## Prediction figures, A 
###########--

prediction.measures <- prediction.columns[!prediction.columns %in% c("hit.selection", "cutoff")]
prediction.labels <- prediction.colnames[!prediction.colnames %in% c("Hit Selection", "Validation Cutoff")]




for (x in 1:length(prediction.measures)){
  prediction.measure.x <- prediction.measures[x]
  prediction.label.x <- prediction.labels[x]
  
  
  
  ## FIgure name
  figure.name <- paste0(prediction.measure.x, "_vs_cutoff")
  figure.title <- paste0(prediction.label.x, " vs. Cuotff")
  
  
  #Set aesthetcis
  hit.selection.colors <- c("black",  "#1A4384")
  hit.selection.colors2 <- c("black",  "firebrick")
  
  
  #Generate figure
  
  prediction.figure <- ggplot(prediction.df, aes(get(prediction.measure.x), cutoff, color=factor(hit.selection)))+
    geom_line(size = 2, alpha = 0.5) +
    geom_point(size = 6)+
    # scale_y_continuous(breaks=c(min(BROAD.tlr.plus.entrez.ordered$ZHL),cutoff,0,max(BROAD.tlr.plus.entrez.ordered$ZHL)),
    #                    labels =round(c(min(BROAD.tlr.plus.entrez.ordered$ZHL),cutoff,0,max(BROAD.tlr.plus.entrez.ordered$ZHL)),1) ) +
    # scale_x_continuous(breaks = c(1, nrow(BROAD.tlr.plus.entrez.ordered) ))+
    scale_color_manual(values = hit.selection.colors2, labels = analysis.labels)+
    labs(title = paste0("CRISPR/Cas9 TNF Validation Screen (Sec Only)    -     ", figure.title), y = "Validation Cutoffs", x = prediction.label.x, color = "Hit Selection Method")+
    # annotate(geom="text", x=0.1*nrow(BROAD.tlr.plus.entrez.ordered), y=0.75*cutoff, label="2.5% cutoff",
    #          color="black", size =7)+
    # annotate(geom="text", x=0.1*nrow(BROAD.tlr.plus.entrez.ordered), y=1.25*cutoff, 
    #          label=paste0(length(which(BROAD.tlr.plus.entrez.ordered$ZHL >= cutoff)), " hits"),
    #          color="black", size =7)+
    theme(axis.text.y = element_text(size=15),
          axis.text.x = element_text(size=15),
          axis.title=element_text(size=20,face="bold"),
          axis.line.y = element_line(colour = 'black', size = 3, lineend = "round"),
          axis.line.x = element_line(colour = 'black', size = 3, lineend = "round"),
          panel.background = element_blank(),   
          panel.grid.major = element_line(colour = "grey"), 
          panel.grid.minor = element_line(colour = "grey", size = 0.25),
          plot.title = element_text(size = 20, hjust = 0.5, vjust = 2, 
                                    margin = margin(b = 6, t = 6)))
  
  setwd(paste0(parnas.wd, "Figures/PredictionFigures/NarrowPrediction_versusCutoff/"))
  ggsave(paste0(figure.name, "_ParnasSecondary_narrowPredict.png"), width = 13.5, height = 12.5)
}




  
  































