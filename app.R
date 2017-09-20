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
Sys.setenv(R_ZIPCMD="/usr/bin/zip")

# global variables
organismAbbr <- NULL
originalHits <- NULL

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
          # selectInput("cutoff_type", "Cutoff type", c("P value", "FDR", "Z score")
          # ),
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
            tabPanel(title = "Enriched Pathways", value = "enrichedPathways", dataTableOutput("enrichedPathways")
            ),
            tabPanel(title = "Gene List", value = "geneList", dataTableOutput("geneList")
            ),
            tabPanel(title = "Network Graph", value = "networkGraph", plotOutput("networkGraph")
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
      # output$cutoffValues <- renderUI({
      #   if (is.null(input$cutoff_type))
      #     return()
      # 
        # switch(input$cutoff_type,
        #        "Z score" = selectInput("cutoff_value", "Cutoff value",
        #                                choices = c("+/-1.96", "+/-2.58"),
        #                                selected = c("+/-1.96")
        #        ),
        #        "P value" = selectInput("cutoff_value", "Cutoff value",
        #                                choices = c("0.05", "0.01"),
        #                                selected = "0.05"
        #        ),
        #        "FDR" = selectInput("cutoff_value", "Cutoff value",
        #                            choices = c("<5%", "<1%"),
        #                            selected = c("<5%")
        #        )
        #   )
      #})
      
      # Read in the input fie  
      output$contents <- renderDataTable({
        inFile <- input$file1
        
        if (is.null(inFile))
          return(NULL)
        
        data <- read.csv(inFile$datapath)
        
        # display the input file dimension
        output$dataSummary <- renderText(c("Your input data have ",  nrow(data),  "rows and ",  ncol(data), "columns."))
   
        return(data)
      })
      
      output$cutoffTypes <- renderUI({
        inFile2 <- input$file1
        if (is.null(inFile2))
          return()
        
        data <- read.csv(inFile2$datapath)
        pulldown_types <- colnames(data)
      
        selectInput("cutoff_type", "Cutoff Types", pulldown_types)
      })
      
      # parameters
      output$organism <- eventReactive(input$goButton, paste("Organism:", input$organism))
      output$pathway  <- eventReactive(input$goButton, paste("Pathway:", input$pathway))
      output$network  <- eventReactive(input$goButton, paste("Network:", input$network))
      output$cutoff_type <- eventReactive(input$goButton, paste("Cutoff Type:", input$cutoff_type))
      output$cutoff_valueH <- eventReactive(input$goButton, paste("High-conf Cutoff Value:", input$cutoff_valueH))
      output$cutoff_valueM <- eventReactive(input$goButton, paste("Med-conf Cutoff Value:", input$cutoff_valueM))
      
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
        dataDir <- "~/data/"

        organism <- input$organism
        organismAbbr <- ifelse(grepl("human", tolower(organism)), 'hsa', 'mmu')
        
        print(organism)
        source("~/TRIAGE/Rscripts/pathway_iteration.R", local = TRUE)
        #source("~/TRIAGE/Rscripts/Network_iteration_V3.R")
        networkType <- "hSTRINGppi.hi"
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
        #pathway.type <- input$pathway
        
        pathwayData <- read.csv(file = paste0(dataDir, "Pathways/", pathway.type, organism, ".csv"))
        
        #seedName <- "IAMinput_HuTNF_manual.csv"
        #seedName <- "TRIAGEinput_HuTNF_CSAfdr_5percCO.csv"
        #setwd(inputDir)
        
        #siRNA.Score <- read.csv(seedName, stringsAsFactors = F)
        siRNA.Score <- read.csv((input$file1)$datapath, stringsAsFactors = F)
        proxyScore <- input$cutoff_type 
        iteration <- 1
        counter <- TRUE
        originalHits <- siRNA.Score$GeneSymbol[siRNA.Score[[proxyScore]] > input$cutoff_valueH]
        
        while (counter == TRUE) {
          
          Hits <- siRNA.Score$EntrezID[siRNA.Score[[proxyScore]] > input$cutoff_valueH]
            
          nonHits <- setdiff(siRNA.Score$EntrezID, Hits)
          
          outPrefix <- paste(pathway.type, iteration, sep = "_")  
          
          # 1) Contraction - [Pathway Analysis]
          siRNA.Score <- ComputeEnrichment(pathwayData, Hits, nonHits, outPrefix, siRNA.Score, iteration)
          siRNA.Score <- data.frame(siRNA.Score, temp = 0, stringsAsFactors = FALSE)
          kName1 <- paste0("KEGG.class.iteration", iteration)
          kName2 <- paste0("KEGG.", iteration)
          names(siRNA.Score)[names(siRNA.Score) == "temp"] <- kName1
          siRNA.Score[[kName1]][siRNA.Score$KEGG == "Yes" & (siRNA.Score[[proxyScore]] > input$cutoff_valueM)] <- 1
          siRNA.Score[[kName1]][siRNA.Score$KEGG != "Yes" & (siRNA.Score[[proxyScore]] > input$cutoff_valueM)] <- 0.5
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
        
        ## 
        
        
        ## Save the results into utput files in the TRIAGEoutputFiles folder
        out <- merge(siRNA.Score, pathDF, by = "GeneSymbol", all = T)
        #setwd(outDir)
        write.csv(out, file = outputFileName, row.names = F)
        final_enriched_pathway_file <- paste0(pathway.type, "_TRIAGE_enrichment_final", ".csv") 
        write.csv(pathEnrich, file = final_enriched_pathway_file, row.names = F)
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
            
            # Used to add text color to GeneSymbols indicate whether they are in the original hit or identified by TRIAGE
            createLink <- function(val1, val2, val3) {
              sprintf('<a href="http://www.genome.jp/kegg-bin/show_pathway?%s0%s" target="_blank" class="btn btn-primary">%s</a>', val1, val2, val3)
            }

            for (i in 1:nrow(pathEnrich))
            {
              # Add hyperlink to KEGG pathway database
              pathwayName <- pathEnrich[i,][1] 
              pathwayID <- pathwayData$PathwayID[match(pathwayName, pathwayData$PathwayName)]
              pathEnrich[i,][1] <- createLink(organismAbbr, pathwayID, pathwayName)
              
              # Add color to indicate whether a gene is in the original hit gene set
              myGeneSet <- pathEnrich[i,7]
              myGenes <- unlist(strsplit(as.character(myGeneSet), ','))
              myOriginalHits <- paste(unlist(originalHits), collapse=', ')
              myGene <- NULL

              for(j in 1:length(myGenes))
              {
                if(grepl(myGenes[j], myOriginalHits)){
                  myGene <- paste(myGene, fontBlue(myGenes[j]), sep = ",")
                }else{
                  myGene <- paste(myGene, fontRed(myGenes[j]), sep = ",")
                }
              }
              myGene <- substring(myGene, 2)
              pathEnrich[i,7] <- myGene
            }
            return(pathEnrich)
            completed4 <- TRUE
          }, escape = FALSE)
        }
      })
      
      # Create the 'Ranked Genes' tab
      completed2 <- TRUE
      output$geneList <- renderUI({
        if(completed2) {
          updateTabsetPanel(session, "inTabset", selected = "geneList")
          
          output$geneList <- renderDataTable({
            # Select only the hit genes in any iteration
            # number of columns to check based on the number of iterations
            numConditions <- NULL
            
            for (i in 1:(iteration - 1))
            {
              ifelse ((i == (iteration - 1)),
                numConditions <- paste(numConditions, "siRNA.Score$\'KEGG.class.iteration", i, "\' == 1.0", sep=""),
                numConditions <- paste(numConditions, "siRNA.Score$\'KEGG.class.iteration", i, "\' == 1.0 | ", sep="")
              )
            }
            
            subSet <- subset(siRNA.Score, eval(parse(text = numConditions)))
            
            # Add a new row of "Total"
            # totalRow <- matrix(c(rep.int(NA,length(subSet))),nrow=1,ncol=length(subSet))
            # colnames(totalRow) <- colnames(subSet)
            # totalRow[1,"GeneSymbol"] <- "Total"
            # totalRow[1,"KEGG.class.iteration1"] <- sum(subSet$'KEGG.class.iteration1'[subSet$'KEGG.class.iteration1' == 1.0])
            # totalRow[1,"KEGG.class.iteration2"] <- sum(subSet$'KEGG.class.iteration2'[subSet$'KEGG.class.iteration2' == 1.0])
            # totalRow[1,"KEGG.class.iteration3"] <- sum(subSet$'KEGG.class.iteration3'[subSet$'KEGG.class.iteration3' == 1.0])
            # totalRow[1,"KEGG.class.iteration4"] <- sum(subSet$'KEGG.class.iteration4'[subSet$'KEGG.class.iteration4' == 1.0])
            # totalRow[1,"KEGG.class.iteration5"] <- sum(subSet$'KEGG.class.iteration5'[subSet$'KEGG.class.iteration5' == 1.0])
            # newSubset <- rbind(totalRow, subSet)
          })
        }
      })
      # 
      # # Create the 'Network Graph' tab
      # completed3 <- TRUE
      # output$networkGraph <- renderUI({
      #   if(completed3) {
      #     updateTabsetPanel(session, "inTabset", selected = "networkGraph")
      #     
      #     output$networkGraph <- renderPlot({networkGraph})
      #   }
      # })
      # 
      #Create the 'Download' tab
      #completed4 <- TRUE
      #observe({
      #if(completed4) {
        updateTabsetPanel(session, "inTabset", selected = "downloads")
        
        output$downloadFiles <- renderUI({
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
    })
  }
    
  #####################
  # Run the application 
  shinyApp(ui = ui, server = server)
}