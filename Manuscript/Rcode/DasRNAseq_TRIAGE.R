####################---
### Preparing TRIAGE input for RNAseq from Das 2018
######################--
## Sam Katz, Oct. 23, 2020


#Load libraries-----
library(dplyr)
library(ggplot2)
library(DESeq2)
library(ggthemes)


#Working directorites -----
rnaseq.wd <- "~/OneDrive\ -\ National\ Institutes\ of\ Health/TRIAGE/DasRNAseq/"
rawdata.wd <- paste0(rnaseq.wd, "GSE103958_RAW/")
deseq2.wd <- paste0(rnaseq.wd, "DESeq_results/")
das.triage.wd <- paste0(rnaseq.wd, "TRIAGEinput/")
analysis.wd <- "~/OneDrive\ -\ National\ Institutes\ of\ Health/Analyis/"
triage.wd <- paste0(analysis.wd, "TRIAGEscripts/")
TRIAGEoutput.dir <- paste0(rnaseq.wd, "TRIAGEoutput/" )
triage.path <- "/Users/sakatz/OneDrive\ -\ National\ Institutes\ of\ Health/Analyis/TRIAGEscripts/"
function.wd <- "~/OneDrive\ -\ National\ Institutes\ of\ Health/TRIAGE/Functions/"
database.wd <- "~/OneDrive\ -\ National\ Institutes\ of\ Health/Databases/"

#### Set experimental types.
cond <- c("Con-BMDM", "IFN-LPS-1-BM", "IFN-LPS-4-BM", "IFN-LPS-12-BM", "IFN-LPS-24-BM")
replicates <- 3

### Count data-----

setwd(rawdata.wd)
for (c in 1:length(cond)){
  for (r in 1:replicates){
    imported.df <-  read.table(paste0(cond[c], "-", r, ".ReadsPerGene.out.tab"), sep = "\t", row.names=1, header = F, stringsAsFactors = F)
    assign(paste0(cond[c], "-", r, "h.raw"), imported.df[-c(1:4),])
  }
}


genomics.range <- rownames(get(paste0(cond[1], "-", 1, "h.raw")))

countdata <- data.frame(matrix(nrow = length(genomics.range), ncol = 0)
                        , row.names = genomics.range
                        )
#Begin id number
id.count.start <- 2786903
id.prefix <- "GSM"

for (c in 1:length(cond)){
  for (r in 1:replicates){
    
    rawcount.df <-  get(paste0(cond[c], "-", r, "h.raw"))
    colnames(rawcount.df) <- paste0(id.prefix, seq.int(id.count.start, id.count.start + 2))
    
    id.count.start <- id.count.start + 3
    
    countdata <- cbind(countdata, rawcount.df)
    
  }
}

#Add 1 to all values
#countdata <- countdata + 1



#### Create matrix input----
setwd(rnaseq.wd)

coldata <- read.csv("lps_ifn_family_matrix.csv", stringsAsFactors = T)
coldata$condition <- NULL

##### dds ----

(dds_test <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~ time))

#### remove technical replicates -----
ddsCollapsed <- collapseReplicates( dds_test,
                                    groupby = dds_test$sample,
                                    run = dds_test$id )

head( as.data.frame( colData(ddsCollapsed)[ ,c("sample","runsCollapsed") ] ), 12 )










######  dds for time subsets ----
timepoints <- c("1hr", "4hr", "12hr", "24hr")
for (t in 1:length(timepoints)){
  
  (dds <- ddsCollapsed[ , ddsCollapsed$time %in% c("0hr", timepoints[t]) ])
  
  dds$time <- droplevels( dds$time )
  
  dds$id <- droplevels( dds$id )
  
  dds$time <- relevel( dds$time, "0hr" )
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds_results <- DESeq(dds)
  
  ##results
  res <- results( dds_results )
  
  #Save the results
  setwd(deseq2.wd)
  write.csv(as.data.frame(res), paste0("LPS_", timepoints[t], "_vs0hr_deseq2.csv"))
  
}


#### Get deseq2 results-----

setwd(deseq2.wd)

for (t in 1:length(timepoints)){
  assign(paste0("deseq2.", timepoints[t], ".df"), read.csv(paste0("LPS_", timepoints[t], "_vs0hr_deseq2.csv"), stringsAsFactors = F))
}


###### figure of data ----
for (t in 1:length(timepoints)){
  FCvPval.plot <- 
    ggplot(get(paste0("deseq2.", timepoints[t], ".df")), aes(y=-1*log2(padj), x=log2FoldChange)) + 
    geom_point()+
    ggtitle(paste0("RNAseq ", timepoints[t], " vs. 0hr: Log FC vs. Log2 adjusted p Value"))+
    #geom_hline(yintercept=-1*log2(0.1), linetype="solid", color = "red")+
    #geom_vline(xintercept= medhits.FC.cutoff, linetype="solid", color = "red")+
    #geom_vline(xintercept= tophits.FC.cutoff, linetype="solid", color = "red")+
    #geom_hline(yintercept = -1*log2(pval.cutoff),  color = "blue")+
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
  
  ggsave(paste0("LPS_", timepoints[t], "_vs0hr_padj_deseq2.png"), width = 12.5, height = 12.5)

  
}

#############highlight cutoffs - padj ------

#Pval cutoff 
padj.cutoff <- 0.05

#Find top hit score and lower score
size_tophits <- 400
size_medhits <- 800

for (t in 1:length(timepoints)){
  screen.df <- get(paste0("deseq2.", timepoints[t], ".df"))
  
  #get dataframe witho only hits above pvalue cutoff
  screen.df.pval <- screen.df[which(screen.df$padj <= 0.05), ]
  
  
  tophits.FC.cutoff <- sort(screen.df.pval$log2FoldChange, decreasing = T)[size_tophits]
  
  medhits.FC.cutoff <- sort(screen.df.pval$log2FoldChange, decreasing = T)[(size_tophits + size_medhits)]
  
  #create figure
  
  FCvPval.cutoff.plot <- 
    ggplot(get(paste0("deseq2.", timepoints[t], ".df")), aes(y=-1*log2(padj), x=log2FoldChange)) + 
    geom_point()+
    ggtitle(paste0("RNAseq ", timepoints[t], " vs. 0hr: Log FC vs. Log2 adjusted p Value"))+
    #geom_hline(yintercept=-1*log2(0.1), linetype="solid", color = "red")+
    geom_vline(xintercept= medhits.FC.cutoff, linetype="solid", color = "red")+
    geom_vline(xintercept= tophits.FC.cutoff, linetype="solid", color = "red")+
    geom_hline(yintercept = -1*log2(padj.cutoff),  color = "blue")+
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
  
  
  setwd(deseq2.wd)
  ggsave(paste0("Cutoffs_", timepoints[t], "_vs0hr_padj_deseq2.png"), width = 12.5, height = 12.5)
  
  
  #Create TRIAGE input
  
  setwd(das.triage.wd)
  
  screen.df.triageinput <- screen.df %>%
    mutate(TRIAGEscore = ifelse(padj <= padj.cutoff & log2FoldChange >= tophits.FC.cutoff, 1, 
                                ifelse(padj <= padj.cutoff & log2FoldChange >= medhits.FC.cutoff, 0.5, 0)))
  
  colnames(screen.df.triageinput)[which(names(screen.df.triageinput) == "X")] <- "GeneSymbol"
  
  
  write.csv(screen.df.triageinput, paste0("TRIAGEinput_", timepoints[t], "_vs0hr_padjCutoff_deseq2.csv"))

}


#############highlight cutoffs - pvalue ------

# #Pval cutoff 
# pval.cutoff <- 0.05
# 
# #Find top hit score and lower score
# size_tophits <- 333
# size_medhits <- 666
# 
# for (t in 1:length(timepoints)){
#   screen.df <- get(paste0("deseq2.", timepoints[t], ".df"))
#   
#   #get dataframe witho only hits above pvalue cutoff
#   screen.df.pval <- screen.df[which(screen.df$pvalue <= 0.05), ]
#   
#   #Round to nearest integer
#   tophits.FC.cutoff <- sort(screen.df.pval$log2FoldChange, decreasing = T)[size_tophits]
#   
#   #round to nearest 0.5 increment
#   medhits.FC.cutoff <- sort(screen.df.pval$log2FoldChange, decreasing = T)[(size_tophits + size_medhits)]
#   
#   #create figure
#   
#   FCvPval.cutoff.plot <- 
#     ggplot(get(paste0("deseq2.", timepoints[t], ".df")), aes(y=-1*log2(pvalue), x=log2FoldChange)) + 
#     geom_point()+
#     ggtitle(paste0("RNAseq ", timepoints[t], " vs. 0hr: Log FC vs. Log2 adjusted p Value"))+
#     #geom_hline(yintercept=-1*log2(0.1), linetype="solid", color = "red")+
#     geom_vline(xintercept= medhits.FC.cutoff, linetype="solid", color = "red")+
#     geom_vline(xintercept= tophits.FC.cutoff, linetype="solid", color = "red")+
#     geom_hline(yintercept = -1*log2(pval.cutoff),  color = "blue")+
#     theme(axis.text.y = element_text(size=20),
#           axis.text.x = element_text(size=20),
#           axis.title=element_text(size=25,face="bold"),
#           axis.line.y = element_line(colour = 'black', size = 0.5),
#           axis.line.x = element_line(colour = 'black', size = 0.5),
#           panel.background = element_blank(),   
#           #panel.grid.major = element_line(colour = "grey"), 
#           #panel.grid.minor = element_line(colour = "grey", size = 0.25),
#           plot.title = element_text(size = 30, hjust = 0.5, vjust = 2, 
#                                     margin = margin(b = 6, t = 6)))
#   
#   
#   setwd(deseq2.wd)
#   ggsave(paste0("Cutoffs_", timepoints[t], "_vs0hr_pvalue_deseq2.png"), width = 12.5, height = 12.5)
#   
#   
#   #Create TRIAGE input
#   
#   setwd(das.triage.wd)
#   
#   screen.df.triageinput <- screen.df %>%
#     mutate(TRIAGEscore = ifelse(pvalue <= pval.cutoff & log2FoldChange >= tophits.FC.cutoff, 1, 
#                                 ifelse(pvalue <= pval.cutoff & log2FoldChange >= medhits.FC.cutoff, 0.5, 0)))
#   
#   colnames(screen.df.triageinput)[which(names(screen.df.triageinput) == "X")] <- "GeneSymbol"
#   
#   
#   write.csv(screen.df.triageinput, paste0("TRIAGEinput_", timepoints[t], "_vs0hr_pvalCutoff_deseq2.csv"))
#   
# }


########## Segmented data - pub figure ------
x.axis.spacer <- 0.1
y.axis.spacer <- 0.5


tophits.color <- "#DF8731"
medhits.color <- "#6C44B2"

tophits.color.light <- "#efc297"
medhits.color.light <- "#ad95d7"

line.size <- 1.5
line.type <- "dashed"





for (t in 1:length(timepoints)){
  screen.df <- get(paste0("deseq2.", timepoints[t], ".df"))
  
  #get dataframe witho only hits above pvalue cutoff
  screen.df.pval <- screen.df[which(screen.df$padj <= 0.05), ]
  
  
  tophits.FC.cutoff <- sort(screen.df.pval$log2FoldChange, decreasing = T)[size_tophits]
  
  medhits.FC.cutoff <- sort(screen.df.pval$log2FoldChange, decreasing = T)[(size_tophits + size_medhits)]
  
  
  
  # ymax.medhits <- -1*log2(min(parnas.df.triageinput[which(parnas.df.triageinput$TRIAGEscore == 0.5), ][["pvalue"]]))
  ymax.tophits <- -1*log2(min(screen.df[which(screen.df$log2FoldChange >= 0), ][["padj"]]))
  
  xmax.tophits <- max(screen.df[which(screen.df$log2FoldChange >= tophits.FC.cutoff), ][["log2FoldChange"]])
  # xmax.medhits <- max(screen.df[which(screen.df$log2FoldChange == 0.5), ][["log2FoldChange"]])
  
  
  y.breaks <- c(0, -1*log2(padj.cutoff), 250, 500, 750, 100)
  y.labels <- as.character(round(y.breaks, 1))
  
  x.breaks <- c(-5, 0, medhits.FC.cutoff, tophits.FC.cutoff, 5, 10)
  x.labels <- as.character(round(x.breaks, 1))
  
  #create figure
  
  FCvPval.cutoff.plot <- 
    ggplot(get(paste0("deseq2.", timepoints[t], ".df"))) + 
    geom_rect(aes(xmin=tophits.FC.cutoff - x.axis.spacer, xmax= xmax.tophits + x.axis.spacer, ymin=(-1*log2(padj.cutoff)), ymax=ymax.tophits + y.axis.spacer), fill = tophits.color.light, color=tophits.color.light, alpha=1, size = 0)+
    geom_rect(aes(xmin=medhits.FC.cutoff, xmax= tophits.FC.cutoff, ymin=(-1*log2(padj.cutoff)), ymax=ymax.tophits + y.axis.spacer), fill = medhits.color.light, color=medhits.color.light, alpha=1, size = 0)+
    geom_point(aes(y=-1*log2(padj), x=log2FoldChange), size = 4)+
    ggtitle(paste0("RNAseq ", timepoints[t], " vs. 0hr: Log FC vs. Log2 adjusted p Value"))+
    geom_hline(yintercept=-1*log2(padj.cutoff), color = "red", linetype = line.type, size = line.size)+
    geom_segment(aes(x = medhits.FC.cutoff, y = -Inf, xend = medhits.FC.cutoff, yend = ymax.tophits + y.axis.spacer), color = "red", linetype = line.type, size = line.size)+
    geom_segment(aes(x = tophits.FC.cutoff, y = -Inf, xend = tophits.FC.cutoff, yend = ymax.tophits + y.axis.spacer), color = "red", linetype = line.type, size = line.size)+
    #geom_hline(yintercept=-1*log2(0.1), linetype="solid", color = "red")+
    # geom_vline(xintercept= medhits.FC.cutoff, linetype="solid", color = "red")+
    # geom_vline(xintercept= tophits.FC.cutoff, linetype="solid", color = "red")+
    # geom_hline(yintercept = -1*log2(pval.cutoff),  color = "blue")+
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
  
  
  
  setwd(deseq2.wd)
  ggsave(paste0("Cutoffs_", timepoints[t], "_vs0hr_pvalue_manuscript.png"), width = 12.5, height = 12.5)
}

########## Segmented data - website figure ------
x.axis.spacer <- 0.1
y.axis.spacer <- 0.5


tophits.color <- "#DF8731"
medhits.color <- "#6C44B2"

tophits.color.light <- "#efc297"
medhits.color.light <- "#ad95d7"

line.size <- 2.5
line.type <- "dashed"





for (t in 1:1){
  screen.df <- get(paste0("deseq2.", timepoints[t], ".df"))
  
  #get dataframe witho only hits above pvalue cutoff
  screen.df.pval <- screen.df[which(screen.df$padj <= 0.05), ]
  
  
  tophits.FC.cutoff <- sort(screen.df.pval$log2FoldChange, decreasing = T)[size_tophits]
  
  medhits.FC.cutoff <- sort(screen.df.pval$log2FoldChange, decreasing = T)[(size_tophits + size_medhits)]
  
  
  
  # ymax.medhits <- -1*log2(min(parnas.df.triageinput[which(parnas.df.triageinput$TRIAGEscore == 0.5), ][["pvalue"]]))
  ymax.tophits <- -1*log2(min(screen.df[which(screen.df$log2FoldChange >= 0), ][["padj"]]))
  
  xmax.tophits <- max(screen.df[which(screen.df$log2FoldChange >= tophits.FC.cutoff), ][["log2FoldChange"]])
  # xmax.medhits <- max(screen.df[which(screen.df$log2FoldChange == 0.5), ][["log2FoldChange"]])
  
  
  y.breaks <- c(0, -1*log2(padj.cutoff), 250, 500, 750, 100)
  y.labels <- as.character(round(y.breaks, 1))
  
  x.breaks <- c(-5, 0, medhits.FC.cutoff, tophits.FC.cutoff, 5, 10)
  x.labels <- as.character(round(x.breaks, 1))
  
  #create figure
  
  FCvPval.cutoff.plot <- 
    ggplot(get(paste0("deseq2.", timepoints[t], ".df"))) + 
    #geom_rect(aes(xmin=tophits.FC.cutoff - x.axis.spacer, xmax= xmax.tophits + x.axis.spacer, ymin=(-1*log2(padj.cutoff)), ymax=ymax.tophits + y.axis.spacer), fill = tophits.color.light, color=tophits.color.light, alpha=1, size = 0)+
    #geom_rect(aes(xmin=medhits.FC.cutoff, xmax= tophits.FC.cutoff, ymin=(-1*log2(padj.cutoff)), ymax=ymax.tophits + y.axis.spacer), fill = medhits.color.light, color=medhits.color.light, alpha=1, size = 0)+
    geom_point(aes(y=-1*log2(padj), x=log2FoldChange), size = 4)+
    ggtitle(paste0("RNAseq ", timepoints[t], " vs. 0hr: Log FC vs. Log2 adjusted p Value"))+
    geom_hline(yintercept=-1*log2(padj.cutoff), color = "red", linetype = line.type, size = line.size)+
    geom_segment(aes(x = medhits.FC.cutoff, y = -Inf, xend = medhits.FC.cutoff, yend = ymax.tophits + y.axis.spacer), color = "red", linetype = line.type, size = line.size)+
    geom_segment(aes(x = tophits.FC.cutoff, y = -Inf, xend = tophits.FC.cutoff, yend = ymax.tophits + y.axis.spacer), color = "red", linetype = line.type, size = line.size)+
    #geom_hline(yintercept=-1*log2(0.1), linetype="solid", color = "red")+
    # geom_vline(xintercept= medhits.FC.cutoff, linetype="solid", color = "red")+
    # geom_vline(xintercept= tophits.FC.cutoff, linetype="solid", color = "red")+
    # geom_hline(yintercept = -1*log2(pval.cutoff),  color = "blue")+
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
  
  
  
  setwd(deseq2.wd)
  ggsave(paste0("Cutoffs_", timepoints[t], "_vs0hr_pvalue_web.png"), width = 12.5, height = 12.5)
}




####### TRIAGE ANALYSIS of time points-----



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


#get input dataframe



for (t in 1:length(timepoints)){
  setwd(das.triage.wd)
  
  triage.input.df <- read.csv(paste0("TRIAGEinput_", timepoints[t], "_vs0hr_padjCutoff_deseq2.csv"), stringsAsFactors = F)
  
  data <- triage.input.df
  
  #Selected_STRINGnetwork.igraph <- G
  message("Networks Loaded")
  
  #message(networkType)
  use.only.commnected.components <- c('Yes')
  
  ## Name of analysis for saved files
  AnalysisTitle <- paste0("RNAseq_", timepoints[t], "_vs0hr_TRIAGE")
  
  print(AnalysisTitle)
  
  
  ## get genome background
  
  library('org.Mm.eg.db')
  x <- org.Mm.egSYMBOL2EG
  backgroundGenes <- read.table(file= "~/TRIAGE/app/data/Mouse_genes_with_protein_product_EntrezID_geneSymbol_lookup.txt", sep="\t", header=TRUE)
  
  
  ### Add EntrezID
  
  
  
  mapped_genes <- mappedkeys(x)
  overlappingGenes <- intersect(as.character(as.list(mapped_genes)), as.character(data$GeneSymbol))
  
  xx <- as.list(x[overlappingGenes])
  y <- unlist(xx)
  y <- data.frame(GeneSymbol = names(y), EntrezID = y, row.names = NULL, stringsAsFactors=FALSE)
  
  tempData <- merge(x=data, y=y, by ="GeneSymbol")
  
  data <- data.frame(tempData, stringsAsFactors = F)
  
  
  #Get genome background
  
  # Get background genes that are not in the input data
  df_backgroundGenes <<- backgroundGenes[!backgroundGenes$GeneSymbol %in% data$GeneSymbol,]
  
  xx <- as.list(x[overlappingGenes])
  y <- unlist(xx)
  y <- data.frame(GeneSymbol = names(y), EntrezID = y, row.names = NULL, stringsAsFactors=FALSE)
  numGeneInInput <- nrow(data)
  
  
  ####################--
  # Define Data Tiers
  ####################--
  
  data <- data %>%
    mutate(ConfidenceCategory = ifelse(TRIAGEscore == 1, "HighConf", ifelse(TRIAGEscore == 0.5, "MedConf", "LowConf")))
  
  df_background <- data.frame(EntrezID = df_backgroundGenes$EntrezID, GeneSymbol=df_backgroundGenes$GeneSymbol, ConfidenceCategory = rep("LowConf", nrow(df_backgroundGenes)))
  # combined input data and background data
  myList <- list(data, df_background)
  data <- data.table::rbindlist(myList, fill = TRUE)
  
  

  
  
  
  ##########################--
  ##### Pre TRIAGE Analysis
  ##########################--
  #Get hits as data frames
  ZscoreHits <- data[data$ConfidenceCategory == "HighConf", c("EntrezID", "GeneSymbol")]
  ZscoreBackground <- data[data$ConfidenceCategory != "HighConf", c("EntrezID", "GeneSymbol")]
  
  #Compute Enrichment
  
  source(paste0(triage.wd, "PathwayAnalysis_Function.R"))
  # pathwayData <- KEGG2017_Human_BP # pathway data is defined b species
  Hits <- ZscoreHits$EntrezID
  nonHits <- ZscoreBackground$EntrezID
  
  file.name <- paste0(AnalysisTitle, "_ZscorePathwaysAnalysis.csv")
  ZscorePath <- paste0(rnaseq.wd, "TRIAGEoutput/")
  Zscore_PathwayAnalyis <- PathwayAnalysis(pathwayData, Hits, nonHits, file.name, ZscorePath)
  
  ########################---
  #### TRIAGE Analysis
  ########################---
  ## Name of analysis for saved files
  
  
  print(AnalysisTitle)
  
  #Set temp working directory for generating enrichment files
  tempWD <- paste0(rnaseq.wd, "TRIAGEoutput/", "temp/")
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
  
  ########################---
  ######## Pathway Output
  #########################--
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
  TRIAGEoutput.dir <- paste0(rnaseq.wd, "TRIAGEoutput/" )
  setwd(TRIAGEoutput.dir)
  
  TRIAGE.cond.output.name <- paste0(AnalysisTitle, "_", "TRIAGE_hits.csv")
  Enrichment.cond.output.name <- paste0(AnalysisTitle, "_", "TRIAGE_enrichment.csv")
  
  write.csv(TRIAGEoutput, file = TRIAGE.cond.output.name)
  write.csv(FinalEnrichment.condensed, file = Enrichment.cond.output.name)
  
  
}


################---
##### Get TRIAGE Hits ----
###############---


condensed.columns <- c("GeneSymbol", "EntrezID", "log2FoldChange", "pvalue", "TRIAGEscore", "ConfidenceCategory", "Pathway", "TRIAGEhit")

for (t in 1:length(timepoints)){
  
  setwd(TRIAGEoutput.dir)
  
  triage.hits <- read.csv(paste0("RNAseq_", timepoints[t], "_vs0hr_TRIAGE", "_", "TRIAGE_hits.csv"), stringsAsFactors = F)
  triage.pathways <- read.csv(paste0("RNAseq_", timepoints[t], "_vs0hr_TRIAGE", "_", "TRIAGE_enrichment.csv"), stringsAsFactors = F)
  topscore.pathways <- read.csv(paste0("RNAseq_", timepoints[t], "_vs0hr_TRIAGE", "_", "ZscorePathwaysAnalysis.csv"), stringsAsFactors = F)

  assign(paste0("RNAseq_", timepoints[t], "_vs0hr_TRIAGE", "_", "hits"), triage.hits[, condensed.columns])
  assign(paste0("RNAseq_", timepoints[t], "_vs0hr_TRIAGE", "_", "pathways"), triage.pathways)
  assign(paste0("RNAseq_", timepoints[t], "_vs0hr_TopScore", "_", "pathways"), topscore.pathways)
  
}




#########--
hit.selections <- c("TopScore", "TRIAGE")

#different figure settings
statistical.tests <- c("pVal", "pValFDR")

#Cutoffs for significance based on test
pValFDR.cutoff <- 0.1
pVal.cutoff <- 0.05

for (s in 1:length(statistical.tests)) {
  statistical.test <- statistical.tests[s]
  
  
  ### Get union of significant pathways
  
  #empty matrix
  sig.pathways <- matrix()
  
  for (t in 1:length(timepoints)){
    
    triage.enrich <- get(paste0("RNAseq_", timepoints[t], "_vs0hr_TRIAGE", "_", "pathways"))
    
    sig.pathways <- c(sig.pathways, triage.enrich$Pathway[which(triage.enrich$pVal <= 0.05)])
    
    sig.pathways <- unique(sig.pathways)
  }
  
  
  #String of metabolic pathways
  metabolic.pathways <- sig.pathways[which(grepl("metabolism", sig.pathways) == TRUE)]
  immune.pathways <- c("Toll-like receptor signaling pathway", "NF-kappa B signaling pathway", "TNF signaling pathway")
  
  for (h in 1:length(hit.selections)){
    hit.selection <- hit.selections[h]
    
    # Assign column with screen name to each data frame
    for(t in 1:length(timepoints)){
      temp.df <- get(paste0("RNAseq_", timepoints[t], "_vs0hr_", hit.selection, "_", "pathways"))
      
      temp.df$ScreenName <-paste0(timepoints[t], "_vs0hr")
      
      temp.df$HitSelection <-paste0(hit.selection)
      
      assign(paste0("dataframe", t), temp.df)
      
    }
    
    #Define statistical test
    
    assign(paste0(hit.selection, "_totalEnrichment"), rbind(get(paste0("dataframe", 1)), 
                                                            get(paste0("dataframe", 2)), 
                                                            get(paste0("dataframe", 3)),
                                                            get(paste0("dataframe", 4)))
    )
  }
  
  
  #remove hit genes names column from top score dataframe
  TopScore_totalEnrichment$HitGeneNames <- NULL
  
  AllScreenEnrichment <- rbind(TopScore_totalEnrichment, TRIAGE_totalEnrichment[, colnames(TopScore_totalEnrichment)])
  
  
  
  
  ###### Shared Pathway figure ----
  
  SharedPathways <- as.matrix(c(immune.pathways, metabolic.pathways))
  
  
  #Filter dataframe for just shared pathways
  AllScreenEnrichment.filter <- AllScreenEnrichment %>%
    filter(Pathway %in% SharedPathways)
  
  
  Pathways.ordered <- as.matrix(SharedPathways) 
  
  Screens.ordered <- paste0(timepoints, "_vs0hr")
  
  
  #reverse order for plot
  Pathways.ordered.figure <- Pathways.ordered[nrow(Pathways.ordered):1,]
  
  
  ###Create figure
  
  
  figure <- ggplot(AllScreenEnrichment.filter, aes(x=factor(ScreenName, levels = Screens.ordered), y=factor(Pathway, levels = Pathways.ordered.figure), size= HitGenes / Genes , fill=get(statistical.test))) +
    geom_point(
      colour="black",
      pch=21,
      alpha = 1
      , stroke = 1
    ) +
    labs(title = paste0("Immuno & Metabolic Pathways: ", statistical.test))+
    #facet_grid(~ScreenName) +
    facet_wrap(~ HitSelection) +
    scale_fill_gradient(low = "#F21A00",  high = "#ba4638", space = "Lab", limit = c(0, get(paste0(statistical.test, ".cutoff"))), na.value = "white")+
    scale_size(range = c(4, 25))+
    theme_light()+
    theme(axis.text.y = element_text(size=5),
          axis.text.x = element_text(size=10, colour = "black", angle = 0, hjust = 0.5),
          axis.title=element_text(size=15,face="bold"),
          axis.text.y.left = element_text(size=5),
          axis.title.x = element_text(vjust = -1),
          aspect.ratio = 2.25/1,
          panel.grid.major = element_line(colour = "grey28"), 
          panel.grid.minor = element_line(colour = "grey28", size = 0.25),
          panel.border = element_rect(fill = NA, 
                                      colour = "black", size = rel(8)))
  
  print(figure)
  setwd(deseq2.wd)
  ggsave(paste0("Pathways_", statistical.test, "_LPSIFNrnaseq.png"), width = 12.5, height = 12.5)
  
  
}

  
  

    
    
    
    
  
  
  
  













  



















































