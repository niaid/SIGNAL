#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(shinyBS)
library(readr)
library(dplyr)
library(stringi)
library(DT)

Sys.setenv(R_ZIPCMD="/usr/bin/zip")

# global variables
organism <- NULL
organismAbbr <- NULL
originalHits <- NULL
networkType <- NULL
completed2 <- FALSE

if (interactive()) {
  
    ##################################################
    # Define UI for application that draws a histogram
    ui <- fluidPage( 
      # style
      theme = "./css/triage.css",
      
      # use shinyjs
      shinyjs::useShinyjs(), br(),
      
      # Application title
      headerPanel(includeHTML("header.html")),
      
      # Sidebar with a slider input for number of bins 
      sidebarLayout(
        
        sidebarPanel(
          # Parameters to be selected
          selectInput(inputId = "organism",
                      label = "Select your organism:",
                      choices = c("Human", "Mouse")
          ),
          selectInput(inputId = "pathway",
                      label = "Select a pathway to enrich:",
                      choices = c("KEGG", "REACTOME", "Gene_Ontology") 
          ),
          selectInput(inputId = "network",
                      label = "Select a PPI database to use:",
                      choices = c("STRING high conf", "STRING medium conf", "BioGRID") 
          ),
          fileInput(inputId= "file1", 
                    label = 'Choose an inputfile to upload'
          ),
          # cutoff values depending the cutoff method chosen
          uiOutput("cutoffTypes"),
          textInput("cutoff_valueH", "High-conf Cutoff Value", placeholder = "High-conf cutoff"),
          bsPopover("cutoff_valueH", "High confidence cutoff value:", "Please enter a value for high confience cutoff, use \"-\" sign for negative value", placement = "bottom", trigger = "hover", options = NULL),
          textInput("cutoff_valueM", "Med-conf Cutoff Value", placeholder = "Med-conf cutoff"),
          bsPopover("cutoff_valueM", "Medium confidence cutoff value:", "Please enter a value for medium confience cutoff, use \"-\" sign for negative value", placement = "bottom", trigger = "hover", options = NULL),
          actionButton("goButton", "Analyze my data !", icon("angle-double-right"), 
                       style="padding:4px; font-size:120%; color: #fff; background-color: rgb(1, 81, 154); border-color: #2e6da4"),
          width = 3
        ),
        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(id = "inTabset",
            tabPanel(title = "Input", value = "contents", 
                     textOutput("dataSummary"), 
                     dataTableOutput("contents"), 
                     textOutput("organism"),
                     textOutput("pathway"),
                     textOutput("network"),
                     textOutput("cutoff_type"),
                     textOutput("cutoff_valueH"),
                     textOutput("cutoff_valueM")
            ),
            tabPanel(title = "Status", value = "status", tags$pre(id = "status")
            ),
            tabPanel(title = "Enriched Pathways", value = "enrichedPathways", 
                     dataTableOutput("enrichedPathways")
            ),
            tabPanel(title = "Gene List", value = "geneList", 
                     dataTableOutput("geneList")
            ),
            tabPanel(title = "Network Graph", value = "myNetworkGraph", 
                     h4('Please select your (1-3) pathways for network graph analysis'), hr(),
                     #textInput("mySelection", label="Your selected pathway IDs:"),
                     #uiOutput("submitGraph"),
                     div(style="display:inline-block",textInput(inputId="mySelection", label="Your selected pathway IDs", value = 0.0)),
                     div(style="display:inline-block",uiOutput("submitGraph")),
                     div(style="display:inline-block",uiOutput("link2Graph")),
                     dataTableOutput("myNetworkGraph")
            ),
            tabPanel(title = "Download", value = "downloads",
                     htmlOutput("downloadFiles"),
                     downloadButton('downloadButton', 'Download all files')
            )
          )    
        )
      ),
      
      # Show a footer using the header style
      headerPanel(includeHTML("footer.html"))
    )
    
    ##################################################
    # Define server logic required to draw a histogram
    server <- function(session, input, output) {
      completed4 <- FALSE

      # Read in the input fie  
      output$contents <- renderDataTable({
        inFile <- input$file1
        
        if (is.null(inFile))
          return(NULL)
        
        data <- read.csv(inFile$datapath)
        
        # display the input file dimension
        #output$dataSummary <- renderText(c("Your input data have ",  nrow(data),  "rows and ",  ncol(data), "columns."))
        data <- datatable(data, rownames = FALSE, options = list(paging=TRUE))
        return(data)
      })
      
      output$cutoffTypes <- renderUI({
        inFile2 <- input$file1
        if (is.null(inFile2))
          return()
        
        data2 <- read.csv(inFile2$datapath)
        pulldown_types <- colnames(data2)
      
        selectInput("cutoff_type", "Cutoff Types", pulldown_types)
      })
      
      # parameters
      # output$organism <- eventReactive(input$goButton, paste("Organism:", input$organism))
      # output$pathway  <- eventReactive(input$goButton, paste("Pathway:", input$pathway))
      # output$network  <- eventReactive(input$goButton, paste("Network:", input$network))
      # output$cutoff_type <- eventReactive(input$goButton, paste("Cutoff Type:", input$cutoff_type))
      # output$cutoff_valueH <- eventReactive(input$goButton, paste("High-conf Cutoff Value:", input$cutoff_valueH))
      # output$cutoff_valueM <- eventReactive(input$goButton, paste("Med-conf Cutoff Value:", input$cutoff_valueM))
      
      # Upon job submission, switch to 'status' tab
      output$status <- eventReactive(input$goButton, 'Analysis started!')

      # Perfrom enrichment
      observeEvent(input$goButton, {
        updateTabsetPanel(session, "inTabset", selected = "status")
        withCallingHandlers({
          shinyjs::html("status", "")
        
        # Global variables
        scriptDir <- "~/TRIAGE/Rscripts/"
        inputDir <- "~/TRIAGE/inputOutputs/TRIAGEinputFiles/"
        outputDir <- "~/TRIAGE/inputOutputs/TRIAGEoutputFiles/"
        dataDir <- "~/TRIAGE/data/"

        organism <- input$organism
        organismAbbr <- ifelse(grepl("human", tolower(organism)), 'hsa', 'mmu')
        
        print(organism)
        source("~/TRIAGE/Rscripts/pathway_iteration.R", local = TRUE)
        
        networkType <- ifelse(grepl("human", tolower(organism)), 'hSTRINGppi.hi', 'mSTRINGhi')
        # # Set the network to be used
        # network <- input$network
        # 
        # if(tolower(organism) == 'human'){
        #   if("hSTRINGhi" %in% network){
        #     networkType <- 'hSTRINGppi.hi'
        #   }
        #   if("hSTRINGmed" %in% network){
        #     networkType <- 'hSTRINGppi.med'
        #   }
        # }
        # else if(tolower(organism) == 'mouse'){
        #   if("mSTRINGhi" %in% network){
        #     networkType <- 'mSTRINGppi.hi'
        #   }
        #   if("mSTRINGmed" %in% network){
        #     networkType <- 'mSTRINGppi.med'
        #   }        
        # }
        
        use.only.commnected.components <- c('Yes')
        
        # Clear any file in the output directory
        do.call(file.remove, list(list.files(outputDir, full.names = TRUE)))
        
        # Set the output file name
        outDir <- "~/TRIAGE/inputOutputs/TRIAGEoutputFiles"
        inputFile <- input$file1
        inputFileName <- inputFile$name
        inputFilePrefix = unlist(strsplit(inputFileName, split='.csv', fixed=TRUE))[1]
        outputFileName <- paste0(inputFilePrefix, "_", networkType, "_TRIGEouput_ALL.csv")
        
        # 1) Seed Pathway Analysis
        setwd(outDir)
        pathway.types <- c("KEGG", "Reactome", "Gene_Ontology")
        pathway.type <- pathway.types[1]
        pathwayData <- read.csv(file = paste0(dataDir, "Pathways/", pathway.type, organism, ".csv"))
        
        # Get input file 
        siRNA.Score <- read.csv((input$file1)$datapath, stringsAsFactors = F)
        proxyScore <- input$cutoff_type 
        iteration <- 1
        counter <- TRUE
        
        # Get a copy of the original list of high-confidence genes
        originalHits <- siRNA.Score$GeneSymbol[siRNA.Score[[proxyScore]] >= input$cutoff_valueH]
        
        # Perform iterative TRIAGE analysis
        while (counter == TRUE) {
          
          Hits <- siRNA.Score$EntrezID[siRNA.Score[[proxyScore]] >= input$cutoff_valueH]
            
          nonHits <- setdiff(siRNA.Score$EntrezID, Hits)
          
          outPrefix <- paste(pathway.type, iteration, sep = "_")  
          
          # 1) Contraction - [Pathway Analysis]
          siRNA.Score <- ComputeEnrichment(pathwayData, Hits, nonHits, outPrefix, siRNA.Score, iteration)
          siRNA.Score <- data.frame(siRNA.Score, temp = 0, stringsAsFactors = FALSE)
          kName1 <- paste0("KEGG.class.iteration", iteration)
          kName2 <- paste0("KEGG.", iteration)
          names(siRNA.Score)[names(siRNA.Score) == "temp"] <- kName1
          siRNA.Score[[kName1]][siRNA.Score$KEGG == "Yes" & (siRNA.Score[[proxyScore]] >= input$cutoff_valueM)] <- 1
          siRNA.Score[[kName1]][siRNA.Score$KEGG != "Yes" & (siRNA.Score[[proxyScore]] >= input$cutoff_valueM)] <- 0.5
          names(siRNA.Score)[names(siRNA.Score) == "KEGG"] <- kName2
          
          hit.Genes <- siRNA.Score$EntrezID[siRNA.Score[[kName1]] == 1]
          myOrignalGenes <- siRNA.Score$GeneSymbol[siRNA.Score[[kName1]] == 1]
          
          # 2) Expansion - [Network Analysis]
          system.time(
            source(paste0(scriptDir, "Network_iteration_V3.R"), local = TRUE)
          )
          siRNA.Score <- data.frame(siRNA.Score, temp1 = "No", temp2 = siRNA.Score[[kName1]], stringsAsFactors = FALSE)
          nName1 <- paste0("Network.", iteration)
          nName2 <- paste0("Network.class.iteration", iteration)
          names(siRNA.Score)[names(siRNA.Score) == "temp1"] <- nName1
          names(siRNA.Score)[names(siRNA.Score) == "temp2"] <- nName2
          siRNA.Score[[nName1]][siRNA.Score$EntrezID %in% gNames2] <- "Yes"
          siRNA.Score[[nName2]][siRNA.Score$EntrezID %in% gNames2 & siRNA.Score[[kName1]] > 0] <- 1
          
          if(iteration != 1 && identical(siRNA.Score[[nName2]], siRNA.Score[[paste0("Network.class.iteration", iteration-1)]])) {
            dupCols <- (ncol(siRNA.Score)-1):ncol(siRNA.Score)
            siRNA.Score <- siRNA.Score[, -dupCols]
            counter <- FALSE
          }
      
          #print(paste("iteration: ", iteration, "\n"))
          message(paste("iteration: ", iteration, "\n"))
          proxyScore <- nName2
          iteration <- iteration + 1
        } # end of while loop
        
        ### Append Enrichment Info --------------------------------------------
        pval_threshold <- 0.05
        #pval_threshold <- input$cutoff_value
        
        iterationNum = iteration - 1
        #enrichFileName <- paste0(outPrefix,".Enrichment_", iterationNum, ".csv")
        enrichFileName <- paste0(outPrefix,".Enrichment_", iterationNum, ".csv")
        pathEnrich <- read.csv(enrichFileName, stringsAsFactors = FALSE)
        
        pathEnrich <- pathEnrich[pathEnrich$pVal < pval_threshold, ]
    
        tempL <- strsplit(pathEnrich$HitGeneNames, split = ", ")
        names(tempL) <- pathEnrich$Pathway
        
        library(reshape2)
        tempL <- lapply(seq(tempL), function(i) {
          m <- melt(tempL[i])
          names(m) <- c("GeneSymbol", paste0("Pathway", i))
          return(m)
        })
        
        tempDF <- Reduce(function(x, y) merge(x, y, by = "GeneSymbol", all = T), tempL)
        
        samp <- tempDF[, 2:ncol(tempDF)]
        pathVector <- sapply(seq(nrow(samp)), function(i) unlist(paste(samp[i, which(!is.na(samp[i, ]))], collapse = " ; ")))
        
        pathDF <- data.frame(GeneSymbol = tempDF$GeneSymbol, Pathway = pathVector, stringsAsFactors = F)
        
        ## Save the results into utput files in the TRIAGEoutputFiles folder
        out <- merge(siRNA.Score, pathDF, by = "GeneSymbol", all = T)
        #setwd(outDir)
        write.csv(out, file = outputFileName, row.names = F)
        final_enriched_pathway_file <- paste0(pathway.type, "_TRIAGE_enrichment_final", ".csv") 
        write.csv(pathEnrich, file = final_enriched_pathway_file, row.names = F)

        # Pass the pathway list to build networkGraph
        #View(pathEnrich)
        sigPathways <<- pathEnrich[,c("Pathway", "Genes", "HitGenes")]
        
        completed <- TRUE
      },
      message = function(m) {
        shinyjs::html(id = "status", html = m$message, add = TRUE)
      })
        
      ## Switch to 'Enriched Pathways' tab adn display partial results
      observe({
        if(completed) {
          updateTabsetPanel(session, "inTabset", selected = "enrichedPathways")
          
          # Display the complete table
          #output$results <- renderDataTable({out[100:115,]})
          # Display the pathway table
          
          output$enrichedPathways <- renderDataTable({
            options = list(autoWidth = TRUE, scrollX = TRUE,
                           columnDefs = list(list(width = '200px', targets = c(4,5))))

            # Used to add hyperlink to KEGG pathway
            fontBlue <- function(val) {
              sprintf('<font color=blue>%s</font>',val)
            }
            
            fontRed <- function(val) {
              sprintf('<font color=red>%s</font>',val)
            }
            
            # Add a hyperlink to KEGG pathway
            createLink <- function(val1, val2, val3) {
              sprintf('<a href="http://www.genome.jp/kegg-bin/show_pathway?%s0%s" target="_blank" class="btn btn-primary">%s</a>', val1, val2, val3)
            }
            
            # Add a hyperlink to KEGG mapper
            link2KEGGmapper <- function(organismAbbr, pathwayID,  myGeneLabels, pathwayName) {
              # Check the data submitted by the form using http post method
              #sprintf('<form target="_blank" enctype="multipart/form-data" method="post" action="http://localhost/cgi-bin/display_form_data.cgi">
              # Create a form for each datatable row
              sprintf('<form target="_blank" enctype="multipart/form-data" method="post" action="http://www.genome.jp/kegg-bin/mcolor_pathway">
                       <input type="hidden" name="map" value="%s0%s">
                       <input type="hidden" name="unclassified" value="%s">
                       <input type="hidden" name="s_sample" value="color">
                       <input type="hidden" name="mode" value="color">
                       <input type="hidden" name="reference" value="white">
                       <input type="submit" style="font-face: \'Comic Sans MS\'; font-size: larger; color: teal; background-color: powderblue; border: 0 none;"value="%s"></form>', organismAbbr, pathwayID, myGeneLabels, pathwayName)
            }
            
            # Used to add text color to GeneSymbols indicate whether they are in the original hit or identified by TRIAGE
            for (i in 1:nrow(pathEnrich))
            {
              # Add hyperlink to KEGG pathway database
              #pathwayName <- pathEnrich[i,][1] 
              #pathwayID <- pathwayData$PathwayID[match(pathwayName, pathwayData$PathwayName)]
              #pathEnrich[i,][1] <- createLink(organismAbbr, pathwayID, pathwayName)
              
              # Add color to indicate whether a gene is in the original hit gene set
              # Also create a hyperlink to KEGG mapper for each pathways
              myGeneSet <- pathEnrich[i,7]
              myGenes <- unlist(strsplit(as.character(myGeneSet), ', '))
              myOriginalHits <- paste(unlist(originalHits), collapse=', ')
              myBlueGene <- NULL
              myBlueGeneID <- NULL
              myBlueGeneLabels <- NULL
              myBlueGeneIDs <- NULL
              myRedGene <- NULL
              myRedGeneID <- NULL
              myRedGeneLabels <- NULL
              myRedGeneIDs <- NULL
              
              
              for(j in 1:length(myGenes))
              {
                if(grepl(myGenes[j], myOriginalHits)){
                  myBlueGene <- paste(myBlueGene, fontBlue(myGenes[j]), sep = ",")
                  myBlueGeneID <- siRNA.Score$EntrezID[which(siRNA.Score$GeneSymbol == myGenes[j])]
                  myBlueGeneLabel <- capture.output(cat(myBlueGeneID, "\t#abebc6,blue\t#abebc6,blue"))
                  myBlueGeneLabels <- paste(myBlueGeneLabels, myBlueGeneLabel, sep="\n")
                  myBlueGeneIDs <- capture.output(cat(myBlueGeneIDs, myBlueGeneLabel))
                  #print(paste(myBlueGeneID, myGenes[j], sep =  " "))
                }else{
                  myRedGene <- paste(myRedGene, fontRed(myGenes[j]), sep = ",")
                  myRedGeneID <- siRNA.Score$EntrezID[which(siRNA.Score$GeneSymbol == myGenes[j])]
                  myRedGeneLabel <- capture.output(cat(myRedGeneID, "\t#ddccff,red\t#ddccff,red\n"))
                  myRedGeneLabels <- paste(myRedGeneLabels, myRedGeneLabel, "\n")
                  myRedGeneIDs <- capture.output(cat(myRedGeneIDs, myRedGeneLabel))
                  #print(paste(myRedGeneID, myGenes[j], sep = " "))
                }
              }
              
              # Add hyperlink to KEGG pathway database
              pathwayName <- pathEnrich[i,][1] 
              pathwayID <- pathwayData$PathwayID[match(pathwayName, pathwayData$PathwayName)]
              mapperHeader <- capture.output(cat("#", organismAbbr,	"CLP/CMP\tBlast_phase\tAll"))
              myGeneLabels <- paste(mapperHeader, stri_replace_all_fixed(myBlueGeneLabels, " ", ""), "\n", stri_replace_all_fixed(myRedGeneLabels, " ", ""), sep = "")
              pathEnrich[i,][1] <- link2KEGGmapper(organismAbbr, pathwayID, myGeneLabels, pathwayName)

              # Display the original hits(BLUE) first, followed by the hits picked up by TRIAGE (RED)
              myGene <- paste(myBlueGene, myRedGene, sep = "")
              myGene <- substring(myGene, 2)
              pathEnrich[i,7] <- myGene
            }
            return(pathEnrich)
            completed2 <- TRUE
          }, escape = FALSE)
        }
      })
      
      # Create the 'Gene List' tab
      output$geneList <- renderUI({
          updateTabsetPanel(session, "inTabset", selected = "geneList")

          output$geneList <- renderDataTable({
            # Select only the hit genes in any iteration
            # number of columns to check based on the number of iterations
            numConditions <- NULL

            for (i in 1:(iteration - 1))
            {
              if (i == (iteration - 1)){
                numConditions <- paste(numConditions, "siRNA.Score$\'KEGG.class.iteration", i, "\' == 1 ", sep="")
              }
              else{
                numConditions <- paste(numConditions, "siRNA.Score$\'KEGG.class.iteration", i, "\' == 1 | ", sep="")
              }
            }

            #print(numConditions)
            subSet <- subset(siRNA.Score, eval(parse(text = numConditions)))

            # Add a new row of "Total" of hit genes in each iteration
            totalRow <- subSet[1,]

            for (k in 1:length(totalRow)) {

              if(grepl("KEGG.class.iteration", colnames(totalRow)[k])) {

                # Get the subset of the dataframe with column value equal to 1
                subSet1 <- subset(subSet, subSet[[colnames(totalRow)[k]]] == 1)
                totalRow[1,colnames(totalRow)[k]] <- sum(subSet1[[colnames(totalRow)[k]]])
              }
              else{
                totalRow[1,k] <- NA
              }
            }
            totalRow[1,"GeneSymbol"] <- "Total"
            newSubset <- rbind(totalRow, subSet)

            # Highlight the 'Total' row using formatStyle()
             dat <- datatable(newSubset, rownames = FALSE, options = list(paging=TRUE)) %>%
              formatStyle('GeneSymbol', target = 'row', backgroundColor = styleEqual(c('Total'), c('orange')))
              return(dat)
          })
          message("completed geneList tab")
      })
      
      # Create the 'Network Graph' tab
      output$myNetworkGraph <- renderUI({
        message("Inside NetworkGraph")
        updateTabsetPanel(session, "inTabset", selected = "myNetworkGraph")
        
        colnames(sigPathways) <- c("Pathways", "TotalGenes", "HitGenes")
        #View(sigPathways)
          
        shinyInput <- function(FUN,id,num,...) {
          inputs <- character(num)
          for (i in seq_len(num)) {
            inputs[i] <- as.character(FUN(paste0(id,i),label=NULL,...))
          }
          inputs
        }
        
        rowSelect <- reactive({
          
          rows=names(input)[grepl(pattern = "srows_",names(input))]
          paste(unlist(lapply(rows,function(i){
            if(input[[i]]==T){
              return(substr(i,gregexpr(pattern = "_",i)[[1]]+1,nchar(i)))
            }
          })))
        })
        
        observe({
          updateTextInput(session, "mySelection", value = rowSelect() ,label = "Your selected pathway IDs:" )
          selectedRows <<- rowSelect()
        })   
        
        output$myNetworkGraph = renderDataTable({
          datatable(cbind(Pick=shinyInput(checkboxInput,"srows_",nrow(sigPathways),value=NULL,width=1), sigPathways[, colnames(sigPathways), drop=FALSE]),
                  options = list(orderClasses = TRUE,
                                 lengthMenu = c(2, 12, 18),
                                 paging = FALSE,
                                 
                                 drawCallback= JS(
                                   'function(settings) {
                                   Shiny.bindAll(this.api().table().node());}')
                                 ),selection='none',escape=F)
        })
        
        output$submitGraph <- renderUI({
          actionButton("submitButton", "Create Network Graph!", icon("angle-double-right"), 
                       style="padding:4px; font-size:120%; color: #fff; background-color: rgb(1, 81, 154); border-color: #2e6da4")
        })  

        observeEvent(input$submitButton, {
          if(length(selectedRows) == 0){
            showNotification("You must select at least ONE pathway!", duration = 5,  type=c("error"))
          }
          if(length(selectedRows) > 3){
            showNotification("You selected more than 3 pathways! Please remove the extra pathways before your re-submission", duration = 10, type=c("error"))
          }
          message(selectedRows)
          
          source("~/TRIAGE/Rscripts/Ranking_plusComments_v2.R", local = TRUE)
          
          Generate_NetworkGraph(selectedRows, organism)
          
        })
        
        observeEvent(input$submitButton, {
          if(length(selectedRows) >= 1 & length(selectedRows) <= 3) {
            output$link2Graph <- renderUI({
              #a(href="file:///Users/songj11/TRIAGE/inputOutputs/TRIAGEoutputFiles/Chimera_STRINGHi_MoTNF.hits.html", target="_blank", "Click here to see the resulting network graph") 
              a("Done! Click here to see the resulting network graph",target="_blank",href="http://localhost/Chimera_STRINGHi_MoTNF.hits.html")
            })
          }
        })
      })
      message("Outside NetworkGraph")
      
      # Create the 'Download' tab
      output$downloadFiles <- renderUI({
        updateTabsetPanel(session, "inTabset", selected = "downloads")
        
        outputFiles <- list.files(path = "~/TRIAGE/inputOutputs/TRIAGEoutputFiles/")
        out <- c("<br><b>All files from your TRIAGE analysis for download:</b><br>")
        
        for(i in 1:length(outputFiles)){
          out <- paste(out, outputFiles[i], sep = "<br>")
        }
        out <- paste0(out, "<br><br>")
        HTML(out)
      })
      
      output$downloadButton <- downloadHandler(
        
        filename = function(){
          paste("TRIAGE_analysis_output","zip",sep=".")
        },
        content = function(filename){
          outputFiles <- list.files(path = "~/TRIAGE/inputOutputs/TRIAGEoutputFiles/")
          zip(zipfile=filename, files = outputFiles)
        },
        contentType = "application/zip"
      )
      
      message("Download content completed")      
    })
  }
    
  #####################
  # Run the application 
  shinyApp(ui = ui, server = server)
}