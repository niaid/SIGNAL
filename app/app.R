# TRIAGE app
# edits May 9, 2018
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
library(data.table)
library(igraph)
library(edgebundleR)
library(shinyAce)
library(mailR)
library(rJava)
library(networkD3)
library(visNetwork)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(reshape2)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(crosstalk)
library(htmltools)
Sys.setenv(R_ZIPCMD="/usr/bin/zip")

# global variables
organism <- NULL
organismAbbr <- NULL
originalHits <- NULL
networkType <- NULL

# Function to validation email addresses
isValidEmail <- function(x) {
  grepl("\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>", as.character(x),
        ignore.case=TRUE)
}

# Set the maximum input file size to 3Mb
options(shiny.maxRequestSize = 3*1024^2)

#if (interactive()) {

    ##################################################
    # Define UI for application that draws a histogram
    ui <- fluidPage(
      
      # Capture user access information
      tags$head(
        tags$title("TRIAGE - Throughput Ranking by Iterative Analysis of Genomic Enrichment"),
        tags$script(src="getIP.js")
      ),

      # style
      theme = "./css/triage.css",

      # use shinyjs
      useShinyjs(), br(),

      # Application title
      headerPanel(includeHTML("header.html")),

      # Sidebar with a slider input for number of bins
      sidebarLayout(

        sidebarPanel(
          # Background color of sidebar panel
          #tags$style(".well {background-color:rgb(1, 81, 154); color: white;}"),
          
          # Global site tag (gtag.js) - Google Analytics 
          tags$head(
            HTML("<!-- Global site tag (gtag.js) - Google Analytics -->"),
            tags$script(src="https://www.googletagmanager.com/gtag/js?id=UA-87121203-23"),
            tags$script(HTML("window.dataLayer = window.dataLayer || [];
                              function gtag(){dataLayer.push(arguments);}
                              gtag('js', new Date());
                    
                              gtag('config', 'UA-87121203-23');"
                        ))
          ),

          # Parameters to be selected
          selectInput(inputId = "organism",
                      label = "Select your organism:",
                      choices = c("Human", "Mouse")
          ),
          selectInput(inputId = "pathway",
                      label = "Select a Database for Enrichment Analysis:",
                      choices = c("KEGG: Biological Processes", "KEGG: Disease Pathways", "KEGG: All Pathways")
          ),
          selectInput(inputId = "network",
                      label = "Select Interactions for Network Analysis:",
                      choices = c("Experimental & Database", "Advanced Options")
          ),
          conditionalPanel(condition = "input.network == 'Advanced Options'",
                           checkboxGroupInput(inputId ="STRING_interaction_sources", 
                                              label = "Select One or More Interaction Sources:",
                                              choices = c(Neighborhood = "neighborhood",
                                                          Fusion = "fusion",
                                                          Cooccurence = "cooccurence",
                                                          Coexpression = "coexpression",
                                                          Experimental = "experimental",
                                                          Database = "database",
                                                          Textmining = "textmining"),
                                              inline = F
                                              , textOutput("txt"))
          ),
          selectInput(inputId = "interaction_confidence_cutoff",
                      label = "Interaction Confidence for Network Analysis:",
                      choices = c("Low (>0.15)" = 150, 
                                  "Medium (>0.4)" = 400,
                                  "High (>0.7)" = 700),
                      selected = 400
          ),
          fileInput(inputId= "file1",
                    label = 'Choose an input file to upload',
                    buttonLabel = "Browse...",
                    # Restrict input file types to .txt and .csv files
                    accept=c("txt/csv", "text/comma-separated-values,text/plain", ".csv")
          ),
          # cutoff values depending the cutoff method chosen
          uiOutput("cutoffTypes"),
          textInput("cutoff_valueH", "High-conf Cutoff Value", placeholder = "High-conf cutoff"),
          bsPopover("cutoff_valueH", "High confidence cutoff value:", "Please enter a value for high confience cutoff, use \"-\" sign for negative value", placement = "bottom", trigger = "hover", options = NULL),
          textInput("cutoff_valueM", "Med-conf Cutoff Value", placeholder = "Med-conf cutoff"),
          bsPopover("cutoff_valueM", "Medium confidence cutoff value:", "Please enter a value for medium confience cutoff, use \"-\" sign for negative value", placement = "bottom", trigger = "hover", options = NULL),
          checkboxInput("includeBackground", "Include background genes"),
          bsPopover("includeBackground", "To include all remaining genes that are not on your input gene list as background", placement = "bottom", trigger = "hover", options = NULL),          actionButton("goButton", "Analyze my data",
                       style="padding:4px; font-size:120%; color: #fff; background-color: rgb(1, 81, 154); border-color: #2e6da4"),
          actionButton("refresh", "Reset", icon("undo"),
                       style="padding:4px; font-size:120%; color: #fff; background-color: rgb(1, 81, 154); border-color: #2e6da4"),
          width = 3
        ),
        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(id = "inTabset",
            tabPanel(title = "Input", value = "contents",
                     htmlOutput("spacer1"),
                     dataTableOutput("contents")
            ),
            tabPanel(title = "Enriched Pathways", value = "enrichedPathways",
                     htmlOutput("spacer2"),
                     dataTableOutput("enrichedPathways")
            ),
            tabPanel(title = "Gene Hits", value = "triageHits",
                     htmlOutput("spacer3"),
                     tabsetPanel(id = 'triageHits',
                                 #TRIAGE Hits table
                                 tabPanel(title = "TRIAGE Gene Hits", value = "triageHits",
                                          dataTableOutput("triageHits")),
                                 tabPanel(title = "Gene Hits By Iteration", value="geneList",
                                          dataTableOutput("geneList")),
                                 tabPanel(title="Graph: Gene Hits By Iteration", value="geneHitsByIteration",
                                          plotOutput("geneHitsByIteration")),
                                 tabPanel(title = "Pathway Enrichments", value = "pathwayEnrich.cond",
                                          dataTableOutput("pathwayEnrich.cond"))
                     )
            ),
            tabPanel(title = "Network", value = "myNetworkGraph",
                     h4('Please select your (1-3) pathways for network graph analysis'), hr(),
                     #textInput("mySelection", label="Your selected pathway IDs:"),
                     #uiOutput("submitGraph"),
                     div(style="display:inline-block",textInput(inputId="mySelection", label="Your selected pathway IDs", value = 0.0)),
                     div(style="display:inline-block",uiOutput("submitGraph")),
                     div(style="display:inline-block",uiOutput("link2Graph")),
                     dataTableOutput("myNetworkGraph")
            ),
            tabPanel(title = "PathNet", value = "graphViews",
                htmlOutput("spacer4"),
                tabsetPanel(id = 'igraphViews',
                      ## Display in igrap
                      tabPanel(title="1st Degree Network", value="graphView1",
                               htmlOutput("graphLegend1"),
                               htmlOutput("graphView1i", width = "100%", height = "700px")
                      ),
                      tabPanel(title="2nd Degree Network", value="graphView2",
                               htmlOutput("graphLegend2"),
                               htmlOutput("graphView2i", width = "100%", height = "700px")
                      ),
                      tabPanel(title = "PathNet Table", value = "PathNetTable",
                               dataTableOutput("PathNetTable"))
                )
            ),
            tabPanel(title = "NetworkD3", value = "networkViews",
                htmlOutput("spacer5"),
                tabsetPanel(id = 'networkD3Views',
                      ## Display in networkD3
                      tabPanel(title="1st Degree D3 Network", value="networkView1",
                               forceNetworkOutput("networkView1D3", width = "100%", height = "700px")),

                      tabPanel(title="2nd Degree D3 Network", value="networkView2",
                              forceNetworkOutput("networkView2D3", width = "100%", height = "700px"))

                      # # Display in visNetwork
                      # tabPanel(title="1st Dimension visNet", value="networkView3",
                      #          visNetworkOutput("networkView3vis", width = "100%", height = "700px")),
                      # 
                      # tabPanel(title="2nd Dimension visNet", value="graphView4",
                      #          visNetworkOutput("networkView4vis", width = "100%", height = "700px"))
                )
            ),
            tabPanel(title = "Download", value = "downloads",
                     htmlOutput("spacer6"),
                     htmlOutput("downloadFiles"),
                     downloadButton('downloadButton', 'Download all files')
            ),
            tabPanel(title = "Help", value = "helpUs",
              htmlOutput("spacer7"),
              tabsetPanel(id = 'helpTab',
                  tabPanel(title = "Contact us", value = "contactUS",
                      uiOutput("contactUS")),
                  tabPanel(title = "Documentation", value = "readMe",
                      uiOutput("documentation")),        
                  tabPanel(title = "Updates", value = "changeLog",
                      uiOutput("changeLog"))                          
              )
            )
          ),
          width = 9
        )
      ),

      # Show a footer using the header style
      headerPanel(includeHTML("footer.html"))
    )

    ##################################################
    # Define server logic required to draw a histogram
    server <- function(session, input, output) {
      
      # Add spacers in the tab panel
      output$spacer1 <- renderUI({
        HTML("<BR><BR>")
      })
      output$spacer2 <- renderUI({
        HTML("<BR><BR>")
      })
      output$spacer3 <- renderUI({
        HTML("<BR>")
      })
      output$spacer4 <- renderUI({
        HTML("<BR>")
      })      
      output$spacer5 <- renderUI({
        HTML("<BR>")
      })  
      output$spacer6 <- renderUI({
        HTML("<BR>")
      })  
      output$spacer7 <- renderUI({
        HTML("<BR>")
      })    
      
      # Global environmental variables
      envs <- Sys.getenv()
      env_names <- names(envs)
      
      ## Set up dataDir
      if('SHINY_SERVER_VERSION' %in% env_names){
        dataDir <- '/srv/shiny-server/data/'
      }else{
        dataDir <- "~/TRIAGE/app/data/"
      }  
      
      # Read in the input fie
      output$contents <- renderDataTable({
        inFile <- input$file1

        if (is.null(inFile))
          return(NULL)

        data <- read.csv(inFile$datapath, stringsAsFactors = FALSE, header=TRUE)
        
        # # Check for duplicated GeneSymbols
        # if(anyDuplicated(data$GeneSymbol)){
        #   showModal(modalDialog(title="User Input Errors", HTML("<h4><font color=red>Duplicated GeneSymbols were found! <br><br>Please remove the duplicates and reload your input file.</font><h4>")))
        # }        
        # 
        # # Check for duplicated GeneSymbols
        # if(anyDuplicated(data$EntrezID)){
        #   showModal(modalDialog(title="User Input Errors", HTML("<h4><font color=red>Duplicated EntrezIDs were found! <br><br>Please remove the duplicates and reload your input file.</font><h4>")))
        # } 
        
        ## Complete list of 19191 protein-encoding genes in human genome
        ## ftp://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/Mus_musculus.gene_info.gz
        humanGenes <- read.table(file=paste0(dataDir, "HGNC_19191_genes_with_protein_product_EntrezID_geneSymbole_lookup.txt"), sep="\t", header=TRUE)
        
        ## Complete list of 23504 genes (mRNAs and ncRNAs) in mouse genome
        ## ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_types/gene_with_protein_product.txt
        mouseGenes <- read.table(file=paste0(dataDir, "Mus_musculus_23504_genes_EntrezID_geneSymbol_lookup.txt"), sep="\t", header=TRUE)
        
        # Check to see if eitehr a 'EntrezID' or a 'GeneSymbol' column is in the input file
        if(!("EntrezID" %in% colnames(data)) && !("GeneSymbol" %in% colnames(data))){
          # Input data do not have EntrezID AND GeneSymbol columns
          showModal(modalDialog(
            title=HTML("<h3><font color=#ff0000>Input file format error!</font></h3>"),
            HTML("Your input file does not contain a required column named 'EntrezID' or 'GeneSymbol'. <br>Please fix your input file and try again!"),
            easyClose = TRUE
          ))
          Sys.sleep(5)
          session$reload()
        }else if(("GeneSymbol" %in% colnames(data)) & !("EntrezID" %in% colnames(data))) {
          # Input data have 'GeneSymbol' column
          # Get EntrezID from GeneSymbol
          if(input$organism == "Human"){
            library('org.Hs.eg.db')
            x <- org.Hs.egSYMBOL2EG
            backgroundGenes <- humanGenes
          } else if(input$organism == "Mouse"){
            library('org.Mm.eg.db') 
            x <- org.Mm.egSYMBOL2EG
            backgroundGenes <- mouseGenes
          }
          
          mapped_genes <- mappedkeys(x)
          overlappingGenes <- intersect(as.character(as.list(mapped_genes)), as.character(data$GeneSymbol))
          
          
          # If no overlapping genes found, catch and handle the error
          if(length(overlappingGenes) == 0){
            showModal(modalDialog(
              title=HTML("<h3><font color=#ff0000>Organism - Gene Set MISMATCH!</font></h3>"),
              HTML("The organism you selected and the organism from which the input data generated do not match. Please select the correct organism or a different input file, then try again"),
              easyClose = TRUE
            ))
          }
          
          # Get background genes that are not in the input data 
          df_backgroundGenes <<- backgroundGenes[!backgroundGenes$GeneSymbol %in% data$GeneSymbol,]
          
          xx <- as.list(x[overlappingGenes])
          y <- unlist(xx)
          y <- data.frame(GeneSymbol = names(y), EntrezID = y, row.names = NULL, stringsAsFactors=FALSE)
          numGeneInInput <- nrow(data)
          
          if(length(overlappingGenes) > 0){
            tempData <- merge(x=data, y=y, by="GeneSymbol")
            data <- tempData
          }
          numGeneWithEntrezID <- nrow(data)
          
          # Display a warning if one or more input genes have no matching EntrezID due obsolete GeneSymbol
          if((numGeneInInput - numGeneWithEntrezID) > 0){
            showModal(modalDialog(title="Warning:", HTML("<h3><font color=red>Only"), numGeneWithEntrezID,HTML("/"),numGeneInInput, HTML("GeneSymbols have mapped EntrezIDs and will be used in this analysis!</font><h3><br>"),
                                  HTML("Either check the organism or update your GeneSymbols to match the official <a href='https://www.genenames.org/cgi-bin/symbol_checker' target=_blank>HGNC</a> symbols if you want to include ALL in this analysis.")))
          }
          
          # Switch/reorder 'EntrezID' to the FIRST column
          data = data[, c(ncol(data), 1:(ncol(data) - 1))]
          rm(tempData,x,y,xx)
          
          message("Input file has a 'GeneSymbol' column!")
        }
        else if(("EntrezID" %in% colnames(data)) & !("GeneSymbol" %in% colnames(data))){
          # Input data have 'EntrezID' column (No 'GeneSymbol' column)
          # Get 'GeneSymbol' from 'EntrezID'
          if(input$organism == "Human"){
            library('org.Hs.eg.db')
            x <- org.Hs.egSYMBOL
            backgroundGenes <- humanGenes
          } else if(input$organism == "Mouse"){
            library('org.Mm.eg.db')
            x <- org.Mm.egSYMBOL
            backgroundGenes <- mouseGenes
          }
          
          #mapped_genes <- as.integer(mappedkeys(x))
          mapped_genes <- mappedkeys(x)
          overlappingGenes <- intersect(mapped_genes, data$EntrezID)
          
          # If no overlapping genes found, catch and handle the error
          if(length(overlappingGenes) == 0){
            showModal(modalDialog(
              title=HTML("<h3><font color=#ff0000>Organism - Gene Set MISMATCH!</font></h3>"),
              HTML("The organism you selected and the organism from which the input data were generated do not match. Please select the correct organism or a different input file, then try again"),
              easyClose = TRUE
            ))
          }
          
          # Get background genes that are not in the input data 
          df_backgroundGenes <<- backgroundGenes[!backgroundGenes$EntrezID %in% data$EntrezID,]
          
          xx <- as.list(x[!is.na(overlappingGenes)])
          y <- unlist(xx)
          y <- data.frame(GeneSymbol = y, EntrezID = names(y), row.names = NULL, stringsAsFactors=FALSE)
          
          if(length(overlappingGenes) > 0){
            # Create a dataframe of the input data with both EntrezID and GeneSymbol
            tempData <- merge(x=data, y=y, by="EntrezID")
            data <- tempData
          }
          
          # Switch/reorder 'EntrezID' to the FIRST column
          data = data[, c(ncol(data), 1:(ncol(data) - 1))]
          
          rm(tempData,x,y,xx)
          message("Input file has a 'EntrezID' column!")
          
          # Having both EntrezID and GeneSymbol 
        }else{ 
          if(input$organism == "Human"){  
            backgroundGenes <- humanGenes
          } else if(input$organism == "Mouse"){
            backgroundGenes <- mouseGenes
          }
          
          # Get background genes that are not in the input data 
          df_backgroundGenes <<- backgroundGenes[!backgroundGenes$EntrezID %in% data$EntrezID,]
          
          # Switch/reorder 'EntrezID' to the FIRST column
          data = data[, c(ncol(data), 1:(ncol(data) - 1))]
          message("Input file has both 'EntrezID' and 'GeneSymbol' columns!")
        }
        
        # Populate GeneSymbolcolumn with EntrezIDs if the corresponding GeneSymbols are not available
        for (i in 1:nrow(data)){
          if(is.na(data$GeneSymbol[i])){
            data$GeneSymbol[i] <- data$EntrezID[i]
          }
        }
        # # Populate GeneSymbolcolumn with EntrezIDs if the corresponding GeneSymbols are not available
        # for (i in 1:nrow(data)){
        #   if(grepl('-', data$GeneSymbol[i])){
        #     data$GeneSymbol[i] <- data$EntrezID[i]
        #   }else if(is.na(data$GeneSymbol[i])){
        #     data$GeneSymbol[i] <- data$EntrezID[i]
        #   }
        # }
        
        # Make sure the EntrezIDs are integers
        data$EntrezID <- as.integer(data$EntrezID)
        data <- data[order(data$EntrezID),]
        
        # # Move the GeneSymbol as the first column
        # siRNA.Score <- data %>%
        #   dplyr::select(GeneSymbol, everything())
        
        # Make a copy of the original input data for later use
        siRNA.Score <<- data
        
        # display the input file dimension
        datatable(data, rownames = FALSE, options = list(paging=TRUE))
      })

      output$cutoffTypes <- renderUI({
        inFile2 <- input$file1
        if (is.null(inFile2))
          return(NULL)

        #data2 <- read.csv(inFile2$datapath)
        # Check validity of the input file
        mtry <- try(read.csv(inFile2$datapath, header = TRUE), 
                    silent = TRUE)
        
        if (class(mtry) != "try-error") {
          data2 <- read.csv(inFile2$datapath, header = TRUE)
        } else {
          showModal(modalDialog(title="User Input Errors:", HTML("<h3><font color=red>Invalid input! Please check your input file.</font><h3>")))
          return(NULL)
        }        
        
        pulldown_types <- c("", colnames(data2))

        selectInput("cutoff_type", "Cutoff Type", pulldown_types)
      })

      # # These user information (userName and userEmail) will collected after the modal is closed/submitted
      # values = reactiveValues(userInfo = "",   ## This is the text that will be displayed
      #                         modal_closed=F)  ## This prevents the values$userInfo output from updating until the modal is closed/submitted
      
      observeEvent(input$refresh, {   
        session$reload()
      })
      
      ## Start perfroming enrichment
      observeEvent(input$goButton, {
        
        ## Check if input file and relevant paremater selected
        ## if not, show an error message
        if(is.null(input$file1)){
          showModal(modalDialog(title="User Input Errors:", HTML("<h3><font color=red>No input file selected!</font><h3>")))
          return(NULL)
        }
        if(input$cutoff_type == ""){
          showModal(modalDialog(title="User Input Errors:", HTML("<h3><font color=red>No 'Cutoff Type' selected!</font><h3>")))
          return(NULL)
        }
        
        # Both cutoff values are required
        if(input$cutoff_valueH == ""){
          showModal(modalDialog(title="User Input Errors:", HTML("<h3><font color=red>Please enter 'High-conf Cutoff Value'!</font><h3>")))
          req(input$cutoff_valueH)
        }else if(input$cutoff_valueM == ""){
          showModal(modalDialog(title="User Input Errors:", HTML("<h3><font color=red>Please enter 'Med-conf Cutoff Value'!</font><h3>")))
          req(input$cutoff_valueM)
        }

        # Once cutoff-type selected and two cutoff values entered
        # Remove duplicate EntrezID rows based on the cutoff values
        if((as.numeric(input$cutoff_valueH) - as.numeric(input$cutoff_valueM)) > 0){
          # sort the dataframe in decending order
          siRNA.Score = siRNA.Score[order(siRNA.Score[,'EntrezID'],-siRNA.Score[,input$cutoff_type]),]
          # remove the duplicated EntrezID row that has a smaller value in cutoff_type column 
          siRNA.Score = siRNA.Score[!duplicated(siRNA.Score$EntrezID),]
        }else{
          # sort the dataframe in acending order
          siRNA.Score = siRNA.Score[order(siRNA.Score[,'EntrezID'], siRNA.Score[,input$cutoff_type]),]
          # remove the duplicated EntrezID row that has a larger value in cutoff_type column 
          siRNA.Score = siRNA.Score[!duplicated(siRNA.Score$EntrezID),]
        }
        
        ## Upon job submission, switch to 'status' tab
        # message("switching to status tab")
        # updateTabsetPanel(session, "inTabset", selected = "status")

        #withCallingHandlers({
        #   shinyjs::html("status", "")

        ## Show progress bar
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        
        progress$set(message = "Performing Enrichment....", value = 0)
        
        ## Set up scriptDir, inputDir, outputDir, wwwDir depending on whether this is used
        # as a standalone tool or on AWS webservice
        if('SHINY_SERVER_VERSION' %in% env_names){
          scriptDir <- '/srv/shiny-server/Rscripts/'
        }else{
          scriptDir <- "~/TRIAGE/app/Rscripts/"
        }

        if('SHINY_SERVER_VERSION' %in% env_names){
          inputDir <- '/srv/shiny-server/inputOutputs/TRIAGEinputFiles/'
        }else{
          inputDir <- "~/TRIAGE/app/inputOutputs/TRIAGEinputFiles/"
        }

        if('SHINY_SERVER_VERSION' %in% env_names){
          outputDir <- '/srv/shiny-server/inputOutputs/TRIAGEoutputFiles/'
        }else{
          outputDir <- "~/TRIAGE/app/inputOutputs/TRIAGEoutputFiles/"
        }

        # To keep a copy of html files for iframe to access
        if('SHINY_SERVER_VERSION' %in% env_names){
          wwwDir <<- '/srv/shiny-server/www/'
        }else{
          wwwDir <<- "~/TRIAGE/app/www/"
        }        
        # Get organism name from user input
        organism <- input$organism
        organismAbbr <- ifelse(grepl("human", tolower(organism)), 'hsa', 'mmu')
        print(organism)

        ## Source other codes depending on whether this is used
        # as a standalone tool or on AWS webservice
        if('SHINY_SERVER_VERSION' %in% env_names){
          source("/srv/shiny-server/Rscripts/pathway_iteration.R", local = TRUE)
        }else{
          source("~/TRIAGE/app/Rscripts/pathway_iteration.R", local = TRUE)
        }

        # Create the network iGraph to be used based on user selection
        #Get ui parameters
        input_netowrk <<- input$network
        network_ConfidenceCutoff <<- as.numeric(input$interaction_confidence_cutoff)
        network_InteractionSources <<- input$STRING_interaction_sources
        
        if(input_netowrk == "Experimental & Database") 
          {
          load(paste0("~/TRIAGE/app/data/Networks/String.", tolower(organism), ".experimental.highConf.igraph.Rdata"))
          load(paste0("~/TRIAGE/app/data/Networks/String.", tolower(organism), ".database.highConf.igraph.Rdata"))
          G <- graph.union(get(paste0("String.", tolower(organism), ".experimental.highConf.igraph")), get(paste0("String.", tolower(organism), ".database.highConf.igraph")))
          { 
            if(network_ConfidenceCutoff <= 400)
              {
              load(paste0("~/TRIAGE/app/data/Networks/String.", tolower(organism), ".experimental.midConf.igraph.Rdata"))
              load(paste0("~/TRIAGE/app/data/Networks/String.", tolower(organism), ".database.midConf.igraph.Rdata"))
              G <<- graph.union(G, get(paste0("String.", tolower(organism), ".experimental.midConf.igraph")), get(paste0("String.", tolower(organism), ".database.midConf.igraph")))
              }
            if(network_ConfidenceCutoff <= 150)
              {
              load(paste0("~/TRIAGE/app/data/Networks/String.", tolower(organism), ".experimental.lowConf.igraph.Rdata"))
              load(paste0("~/TRIAGE/app/data/Networks/String.", tolower(organism), ".database.lowConf.igraph.Rdata"))
              G <- graph.union(G, get(paste0("String.", tolower(organism), ".experimental.lowConf.igraph")), get(paste0("String.", tolower(organism), ".database.lowConf.igraph")))
              }
            }
          } else if(input_netowrk == "Advanced Options")
          {
            for (i in 1:length(network_InteractionSources)) {
              load(paste0("~/TRIAGE/app/data/Networks/String.", tolower(organism), ".", network_InteractionSources[i], ".highConf.igraph.Rdata"))
              if(exists("G")) {G <- graph.union(G, get(paste0("String.", tolower(organism), ".", network_InteractionSources[i], ".highConf.igraph")))}
              else {G <- get(paste0("String.", tolower(organism), ".", network_InteractionSources[i], ".highConf.igraph"))}
              }
            {
              if(network_ConfidenceCutoff <= 400)
                {
                for (i in 1:length(network_InteractionSources)) {
                  load(paste0("~/TRIAGE/app/data/Networks/String.", tolower(organism), ".", network_InteractionSources[i], ".midConf.igraph.Rdata"))
                  G <- graph.union(G, get(paste0("String.", tolower(organism), ".", network_InteractionSources[i], ".midConf.igraph")))
                }
                }
              if(network_ConfidenceCutoff <= 150)
                {
                for (i in 1:length(network_InteractionSources)) {
                  load(paste0("~/TRIAGE/app/data/Networks/String.", tolower(organism), ".", network_InteractionSources[i], ".midConf.igraph.Rdata"))
                  G <- graph.union(G, get(paste0("String.", tolower(organism), ".", network_InteractionSources[i], ".midConf.igraph")))
                }
              }
            }
          }
        ##Push graph to global environment
        #G <<- G
        
        #Selected_STRINGnetwork.igraph <- G
        message("Networks Loaded")

        #message(networkType)
        use.only.commnected.components <- c('Yes')

        # Generate logfiles
        # user_login.log - with user input information
        # logDir <- getwd()
        # message(logDir, "***")
        # userLoginInfo <- c(input$userName, input$userEmail)
        # write.table(userLoginInfo, file = paste0(logDir, '/user_login.log'), append = TRUE)

        # user_access.log
        # Capture user access information
        # IP <- reactive({ input$getIP })
        # observe({
        #   cat(capture.output(str(IP()), split=TRUE))
        #   userAccessInfo <- capture.output(str(IP()), split=TRUE)
        #   write.table(userAccessInfo, file = paste0(logDir, '/user_access.log'), append = TRUE, quote = TRUE, sep = " ",
        #               eol = "\n", na = "NA", dec = ".", row.names = TRUE,
        #               col.names = TRUE, qmethod = c("escape", "double"))
        # })

        # Create uesr-specific directory
        # userDir <- input$userEmail
        # userDir <- gsub("\\.", "_", userDir)
        # userDir <- gsub("@", "_", userDir)
        # setwd(outputDir)
        # outDir <- c(getwd(), "/", userDir)

        # Remove all files in the outputDir
        unlink(outputDir, recursive = FALSE)

        # Create user-specific directory
        # if(!dir.exists(userDir)[1]){
        #   dir.create(userDir)
        # }
        # else{
        #   unlink(outDir, recursive = FALSE)
        # }
        # setwd(userDir)
        
        # Create user-specific directory using system time
        userDir <- format(Sys.time(),"%Y%m%d%H%M%S%ms")
        outDir <- paste0(outputDir,"/", userDir)
        dir.create(outDir)
        setwd(outDir)
        
        # Set the output file name
        inputFile <- input$file1
        inputFileName <- inputFile$name
        #inputFilePrefix = (unlist(strsplit(inputFileName, split='.csv', fixed=TRUE)))[1]
        inputFilePrefix <- tools::file_path_sans_ext(inputFileName)

        outputFileName <- paste0(inputFilePrefix, "_", "_TRIAGEouput_ALL.csv")

        # 1) Seed Pathway Analysis
        # if('SHINY_SERVER' %in% env_names){
        #   setwd("/srv/shiny-server/inputOutputs/TRIAGEoutputFiles")
        # }
        # else{
        #   setwd("~/TRIAGE/app/inputOutputs/TRIAGEoutputFiles")
        # }

        # Pathway types
        #pathway.types <- c("KEGG", "Reactome", "Gene_Ontology")
        #pathway.type <- pathway.types[1]
        pathway.type <- input$pathway
        
        #Get appropiate pathway document suffix
        if (pathway.type == 'KEGG: Biological Processes'){
          pathway.type <- 'BiologicalProcesses'
        } 
        if (pathway.type == 'KEGG: Disease Pathways'){
          pathway.type <- 'Disease'
        } 
        if (pathway.type == 'KEGG: All Pathways'){
          pathway.type <- 'All'
        } 
        
        pathwayData <- read.csv(file = paste0(dataDir, "Pathways/KEGG2017_", organism, "_", pathway.type, ".csv"))
        # Get input file
        # siRNA.Score <- read.csv((input$file1)$datapath, stringsAsFactors = F)
        # Populate GeneSymbolcolumn with EntrezIDs if the corresponding GeneSymbols are not available
        # for (i in 1:nrow(siRNA.Score )){
        #   if(siRNA.Score$GeneSymbol[i] == '-'){
        #     siRNA.Score$GeneSymbol[i] <- siRNA.Score$EntrezID[i]
        #   }
        # }        
        
        proxyScore <- input$cutoff_type
        
        ## To include background genes 
        includeBackground <- input$includeBackground
        if(includeBackground){
          # To add a value to cutoffType and the value should be based on values of the two cutoff_values
          if((as.numeric(input$cutoff_valueH) - as.numeric(input$cutoff_valueM)) > 0){
            backgroundValue = as.numeric(input$cutoff_valueM) -  (as.numeric(input$cutoff_valueM)/10)
            df_background <- data.frame(EntrezID = df_backgroundGenes$EntrezID, GeneSymbol=df_backgroundGenes$GeneSymbol, proxyScore = rep(backgroundValue, nrow(df_backgroundGenes)))
          }else{
            backgroundValue = as.numeric(input$cutoff_valueM) +  (as.numeric(input$cutoff_valueM)/10)
            df_background <- data.frame(EntrezID = df_backgroundGenes$EntrezID, GeneSymbol=df_backgroundGenes$GeneSymbol, proxyScore = rep(backgroundValue, nrow(df_backgroundGenes)))
          }
          
          # change the name the cutoff_type from 'proxScore' to what it stands in the original data
          names(df_background)[names(df_background) == 'proxyScore'] <- input$cutoff_type
          
          # combined input data and background data
          myList <- list(siRNA.Score, df_background)
          siRNA.Score <- rbindlist(myList, fill = TRUE)
        }

        cutoffType <<- input$cutoff_type
        cutoffHigh <<- as.numeric(input$cutoff_valueH)
        cutoffMed <<- as.numeric(input$cutoff_valueM)
        iteration <- 1
        counter <- TRUE

        # Get a copy of the original list of high-confidence genes
        originalHits <- siRNA.Score$GeneSymbol[siRNA.Score[[proxyScore]] >= input$cutoff_valueH]

        # Check to see if the genes from the input file match the genes from the pathwayData
        # No match in case the input genes are from Mouse, but the genes from the pathwayData are from Human
        matchingHits <- intersect(siRNA.Score$EntrezID,pathwayData$EntrezID)

        # If zero length, catch and handle the error
        if(length(matchingHits) == 0){
          message("Organism - gene set mismatch! You selected either the incorrect organism or a wrong data file.")
          showNotification("organism-gene set MISMATCH! \nThe organism you selected and the organism from which the input data generated do not match.", duration = NULL,
                           action = a(href = "javascript:location.reload();", "Reload page")
          )
        }

        # Perform iterative TRIAGE analysis
        while (counter == TRUE) {

          ## Show progress
          # Increment the progress bar, and update the detail text.
          progress$inc(1/(iteration*3))
          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(0.1)
          
          Hits <- siRNA.Score$EntrezID[siRNA.Score[[proxyScore]] >= input$cutoff_valueH]

          nonHits <- setdiff(siRNA.Score$EntrezID, Hits)

          outPrefix <- paste("KEGG", iteration, sep = "_")

          # 1) Contraction - [Pathway Analysis]
          message(getwd())
          siRNA.Score <- ComputeEnrichment(pathwayData, Hits, nonHits, outPrefix, siRNA.Score, iteration)
          siRNA.Score <- data.frame(siRNA.Score, temp = 0, stringsAsFactors = FALSE)
          kName1 <- paste0("KEGG.class.iteration", iteration)
          kName2 <- paste0("KEGG.", iteration)
          names(siRNA.Score)[names(siRNA.Score) == "temp"] <- kName1
          siRNA.Score[[kName1]][siRNA.Score$KEGG == "Yes" & (siRNA.Score[[proxyScore]] >= as.numeric(input$cutoff_valueM))] <- as.numeric(input$cutoff_valueH)
          siRNA.Score[[kName1]][siRNA.Score$KEGG != "Yes" & (siRNA.Score[[proxyScore]] >= as.numeric(input$cutoff_valueM))] <- as.numeric(input$cutoff_valueM)
          names(siRNA.Score)[names(siRNA.Score) == "KEGG"] <- kName2

          hit.Genes <- siRNA.Score$EntrezID[siRNA.Score[[kName1]] == input$cutoff_valueH]
          myOrignalGenes <- siRNA.Score$GeneSymbol[siRNA.Score[[kName1]] == input$cutoff_valueH]

          # 2) Expansion - [Network Analysis]
    message("*", paste0(scriptDir, "Network_iteration_V3.R"), "**")
    
          source(paste0(scriptDir, "Network_iteration_V3.R"), local = TRUE)
          siRNA.Score <- data.frame(siRNA.Score, temp1 = "No", temp2 = siRNA.Score[[kName1]], stringsAsFactors = FALSE)
          nName1 <- paste0("Network.", iteration)
          nName2 <- paste0("Network.class.iteration", iteration)
          names(siRNA.Score)[names(siRNA.Score) == "temp1"] <- nName1
          names(siRNA.Score)[names(siRNA.Score) == "temp2"] <- nName2
          siRNA.Score[[nName1]][siRNA.Score$EntrezID %in% gNames2] <- "Yes"
          # siRNA.Score[[nName2]][siRNA.Score$EntrezID %in% gNames2 & siRNA.Score[[kName1]] > 0] <- 1
          siRNA.Score[[nName2]][siRNA.Score$EntrezID %in% gNames2 & siRNA.Score[[kName1]] >= input$cutoff_valueM] <- input$cutoff_valueH

          
          
          
          #if(iteration != 1 && identical(siRNA.Score[[nName2]], siRNA.Score[[paste0("Network.class.iteration", iteration-1)]])) {
          if((iteration != 1 && identical(siRNA.Score[[nName2]], siRNA.Score[[paste0("Network.class.iteration", iteration-1)]])) 
             || (iteration >= 5 
                 && identical(siRNA.Score[[nName2]], siRNA.Score[[paste0("Network.class.iteration", iteration-2)]]) 
                 && identical(siRNA.Score[[paste0("Network.class.iteration", iteration-1)]], siRNA.Score[[paste0("Network.class.iteration", iteration-3)]])
                 && (length(siRNA.Score$EntrezID[siRNA.Score[[nName2]]== 1]) > length(siRNA.Score$EntrezID[siRNA.Score[[paste0("Network.class.iteration", iteration-1)]] == 1])))) { 
           #dupCols <- (ncol(siRNA.Score)-1):ncol(siRNA.Score)      #This just shaves off the last two columns 
            #siRNA.Score <- siRNA.Score[, -dupCols]
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
        pathVector <- sapply(seq(nrow(samp)), function(i) unlist(paste(samp[i, which(!is.na(samp[i, ]))], collapse = ", "))) #Switched to comma seperator from " ; " 

        pathDF <- data.frame(GeneSymbol = tempDF$GeneSymbol, Pathway = pathVector, stringsAsFactors = F)

        ## Save the results into output files in the TRIAGEoutputFiles folder
        out <- merge(siRNA.Score, pathDF, by = "GeneSymbol", all = T)

        message(getwd(), "#####")

        write.csv(out, file = outputFileName, row.names = F)
        final_enriched_pathway_file <- paste0(pathway.type, "_TRIAGE_enrichment_final", ".csv")
        write.csv(pathEnrich, file = final_enriched_pathway_file, row.names = F)

        # Pass the pathway list to build networkGraph
        sigPathways <<- pathEnrich[,c("Pathway", "Genes", "HitGenes")]

        completed <- TRUE
      #},
      # message = function(m) {
      #   shinyjs::html(id = "status", html = m$message, add = TRUE)
      # })
        
        
      ###################################
      ### Generate Condensed Output Files
      ###################################
      
      TRIAGEoutput <- read.csv(outputFileName, stringsAsFactors = F)
        cutoffType <<- input$cutoff_type
        cutoffHigh <<- as.numeric(input$cutoff_valueH)
        cutoffMed <<- as.numeric(input$cutoff_valueM)
        
        #################
        # Set High and Low confidence Label and TRIAGE Label
        ##################
        
        #Set confidence catgeor
        TRIAGEoutput <- TRIAGEoutput %>%
          mutate(ConfidenceCategory = ifelse(get(cutoffType, as.environment(TRIAGEoutput)) >= cutoffHigh, "HighConf",
                                             ifelse(get(cutoffType, as.environment(TRIAGEoutput)) >= cutoffMed & get(cutoffType, as.environment(TRIAGEoutput)) < cutoffHigh,"MedConf",
                                                    "")))
        
        FinalIterationNetworkColumn <- paste0("Network.class.iteration", iterationNum)
        
        TRIAGEoutput <- TRIAGEoutput %>%
          mutate(TRIAGEhit = ifelse(get(FinalIterationNetworkColumn, envir = as.environment(TRIAGEoutput)) == cutoffHigh, 
                                    "Yes",
                                    ""))
        ################
        # Generate Matrices for TRIAGE Hits, high conf hits, med conf hits
        #Get filtered TRIAGEhits                                           # This is where I put together what is considered a "hit" by TRIAGE (IAM). Any gene that had a score of 1 in the last network analysis step
        
        TRIAGEhits <- filter(TRIAGEoutput, TRIAGEhit == "Yes")
        TRIAGEhits.matrix <- matrix(TRIAGEhits$EntrezID)                      # Created a matrix of all the genes that are "hits" 
        
        # Get filtered TRIAGE High/Med Conf
        TRIAGEhits.highConf <- filter(TRIAGEoutput, TRIAGEhit == "Yes" 
                                      & ConfidenceCategory == "HighConf")
        TRIAGEhits.highConf.matrix <- matrix(TRIAGEhits.highConf$EntrezID)
        TRIAGEhits.highConf.matrix.GS <- matrix(TRIAGEhits.highConf$GeneSymbol)
        
        TRIAGEhits.medConf <- filter(TRIAGEoutput, TRIAGEhit == "Yes" 
                                     & ConfidenceCategory == "MedConf")
        TRIAGEhits.medConf.matrix <- matrix(TRIAGEhits.medConf$EntrezID)
        TRIAGEhits.medConf.matrix.GS <- matrix(TRIAGEhits.medConf$GeneSymbol)
        
        #numTotal variabale to be used in GeneList tab to generate graph
        
        numTotal  <<- length(TRIAGEhits.highConf.matrix) + length(TRIAGEhits.medConf.matrix)
        
        #############################################################     # Using the methods from CARD (only swithced to EntrezID over GeneSymbol)
        #            PageRank Algorithm to All Genes
        #############################################################               
        # AVERAGE DUPLICATED ROWS
        siRNA.Score <- TRIAGEhits[which(!is.na(TRIAGEhits$EntrezID)),]    
        OverallDegree <- degree(G)
        Screen_Genes.for.network.analysis <- intersect(siRNA.Score$EntrezID,V(G)$name)
        Graph <- induced.subgraph(G,Screen_Genes.for.network.analysis)
        
        #############################################################
        #            Select sub-graphs from hit genes
        #############################################################
        Subset_Genes.for.network.analysis <- Screen_Genes.for.network.analysis
        
        SubGraph <- induced.subgraph(Graph,Subset_Genes.for.network.analysis)
        
        # if(tolower(use.only.commnected.components) == "yes"){
        #   SubGraph <- induced.subgraph(SubGraph,names(which(degree(SubGraph) > 0)))
        # }
        
        #############################################################
        #            Calculation of Node properties
        #############################################################
        Temp.Articulation <- V(SubGraph)$name[articulation.points(SubGraph)]
        Articulation <- rep(0, length(V(SubGraph)$name))
        Articulation[match(Temp.Articulation, (V(SubGraph)$name))] <- 1
        
        siRNA.Score.Formatted <- siRNA.Score
        
        gNames <- V(SubGraph)$name
        ###############################################################################
        #                 Create Network for D3 Rendering for Hit Genes
        ###############################################################################
        GraphEdgesHitNames <- get.data.frame(SubGraph, what = "edges")
        GraphEdgesHitNames <- GraphEdgesHitNames[!duplicated(GraphEdgesHitNames), ]
        source <- target <- rep(NA, nrow(GraphEdgesHitNames))
        tempGenes <- union(GraphEdgesHitNames$from,GraphEdgesHitNames$to)
        for(i in 1:length(tempGenes)){
          source[which(GraphEdgesHitNames$from == tempGenes[i])] <- i-1
          target[which(GraphEdgesHitNames$to == tempGenes[i])] <- i-1
        }
        
        GraphEdgesHitNumber <- data.frame(source,target)
        
        
        GraphNodesHit <- data.frame(GeneMappingID = rep(0:(length(tempGenes)-1)), EntrezID = tempGenes)
        
        
        ###############################################################################
        #                 Add back GeneSymbol and Hit Designation to output
        ###############################################################################
        GraphNodesHit <- merge(GraphNodesHit, TRIAGEhits[, c("EntrezID", "GeneSymbol", "TRIAGEhit")], 
                               by.x = "EntrezID", by.y = "EntrezID", all.x = T)
        
        #############**************************************************################    # Now getting a data frame for the edges and a data frame for the nodes
        
        ############# EDGE INFO and NODE INFO are sent to GLOBAL ENVIRONMENT to be used by RANKING
        
        EdgeInfo <<- GraphEdgesHitNumber
        NodeInfo <<- GraphNodesHit
        
        ######################################
        
        #Add GeneSymbols to Edge dataframe                                              
        Edge.source <- merge(EdgeInfo, NodeInfo[, c("GeneMappingID", "GeneSymbol")], by.x = "source", by.y = "GeneMappingID", all.x = TRUE)
        
        names(Edge.source)[names(Edge.source)=="GeneSymbol"] <- "source.ID"
        
        Edge.target <- merge(Edge.source, NodeInfo[, c("GeneMappingID", "GeneSymbol")], by.x = "target", by.y = "GeneMappingID", all.x = TRUE)
        
        names(Edge.target)[names(Edge.target)=="GeneSymbol"] <- "target.ID"
        
        #####
        #Merge TRIAGEhits to NodeInfo
        Scores_and_nodes <- merge(NodeInfo[, c("GeneMappingID", "GeneSymbol")],                      #Pairing up the "Node Info" with the gene info (such as gene symbol and groupings)
                                  TRIAGEhits, 
                                  by.x = "GeneSymbol", by.y = "GeneSymbol", all.y = T)
        ###################
        ######Create Data frame with Pathway interactome
        #Convert Target values to pathways
        EdgeInfo.TargetPathways <- merge(EdgeInfo, Scores_and_nodes[, c("GeneMappingID", "Pathway")],
                                         by.x = "target", by.y = "GeneMappingID", all.x = T)
        #Remove NAs
        EdgeInfo.TargetPathways <- na.omit(EdgeInfo.TargetPathways)
        
        #Aggregate
        EdgeInfo.TargetPathways_Sum <- EdgeInfo.TargetPathways %>% 
          group_by(source) %>% summarise(Pathway = toString(Pathway))
        
        
        
        #Convert Source values to pathways
        EdgeInfo.SourcePathways <- merge(EdgeInfo, Scores_and_nodes[, c("GeneMappingID", "Pathway")],
                                         by.x = "source", by.y = "GeneMappingID", all.x = T)
        #Remove NAs
        EdgeInfo.SourcePathways <- na.omit(EdgeInfo.SourcePathways)
        
        #Aggregate
        EdgeInfo.SourcePathways_Sum <- EdgeInfo.SourcePathways %>% 
          group_by(target) %>% summarise(Pathway = toString(Pathway))
        
        
        
        ########## Combine data frames
        #Align column names and stack data frames                                                         #The analysis assumed directionality of interactions, but we're ignoring it here, so combining the "target" and Source" to one dataframe.
        colnames(EdgeInfo.SourcePathways_Sum) <- c("source", "Pathway")
        EdgePathways_stacked <- rbind(EdgeInfo.TargetPathways_Sum, EdgeInfo.SourcePathways_Sum)
        
        
        #Aggregate
        EdgePathways_stacked_Sum <- EdgePathways_stacked %>% 
          group_by(source) %>% summarise(Pathway = toString(Pathway))
        
        ##Add counts to pathways (number of genes in list that are part of each pathway)
        EdgePathways_stacked_Sum$Pathway.counts <- NA
        
        for (i in 1:length(EdgePathways_stacked_Sum$Pathway)) {
          temp.string <- unlist(strsplit(EdgePathways_stacked_Sum$Pathway[i], ", "))
          for (j in 1:length(temp.string)) {
            name.j <- temp.string[j]
            count.j <- length(grep(name.j, temp.string, fixed = T))
            out <- paste0(name.j, " (", count.j, ")")
            EdgePathways_stacked_Sum$Pathway.counts[i] <-  ifelse(j == 1, out, paste(out, EdgePathways_stacked_Sum$Pathway.counts[i], sep = ", ")) 
          }
        }
        
        ########Remove Duplicates of Pathway names
        
        EdgePathways_stacked_Sum$NetworkGenePathways <- sapply(strsplit(EdgePathways_stacked_Sum$Pathway.counts, ", ", fixed = TRUE), function(x) 
          paste(unique(x), collapse = ", "))
        
        ###########Order from pathway with highest number of genes to lowest
        for (k in 1:length(EdgePathways_stacked_Sum$Pathway)) {
          pathway.string <- unlist(strsplit(EdgePathways_stacked_Sum$NetworkGenePathways[k], ", "))
          pathway.values <- as.numeric(gsub("[\\(\\)]", "", regmatches(pathway.string, gregexpr("\\(.*?\\)", pathway.string))))
          names(pathway.values) <- order(pathway.values, decreasing = T)
          EdgePathways_stacked_Sum$NetworkGenePathways[k] <- paste(pathway.string[as.numeric(names(pathway.values))], collapse = ", ")
        }
        
        
        #######Add in entrezID for source genemappings
        EdgePathways_stacked_EntrezID <- merge(Scores_and_nodes[, c("GeneMappingID", "EntrezID")], EdgePathways_stacked_Sum,
                                               by.x = "GeneMappingID", by.y = "source",
                                               all.y = T)
        
        #######Add to triage Out put file
        TRIAGEoutput <- merge(TRIAGEoutput, EdgePathways_stacked_EntrezID[, c("EntrezID", "NetworkGenePathways")],
                              by.x = "EntrezID", by.y = "EntrezID",
                              all.x = T)
        
        ##########################
        # Add interacting genes to output files
        
        
        #Aggregate the edges to be summarised to each genemap ID                                            #pulling together all genes that interactect with a specifc gene and putting them in one row seprated by comma
        Edge_source_summary <- aggregate(target.ID ~ source.ID, data = Edge.target, paste, collapse = ", ")
        Edge_target_summary <- aggregate(source.ID ~ target.ID, data = Edge.target, paste, collapse = ", ")
        
        
        #Align column names and stack data frames                                                         #The analysis assumed directionality of interactions, but we're ignoring it here, so combining the "target" and Source" to one dataframe.
        colnames(Edge_target_summary) <- c("source.ID", "target.ID")
        Edge_summary_stacked <- rbind(Edge_source_summary, Edge_target_summary)
        
        #Aggregate stacked data to get unique values
        Edge_summary <- aggregate(target.ID ~ source.ID, data = Edge_summary_stacked, paste, collapse = ", ")
        
        #Update Names                                                                                     #Now a dataframe is being created where each gene selected by TRIAGE (or part of the highlighted groups) has a list "Ntwrk.all" that lists all other genes from TRAIGE that it is predicted to interact with.
        colnames(Edge_summary) <- c("GeneSymbol", "InteractingGenes")
        
        #Merge with scores 
        TRIAGEoutput <- merge(TRIAGEoutput, Edge_summary, by.x = "GeneSymbol", by.y = "GeneSymbol", all.x = T)
        
        ####################
        ### Create Condensed Output File
        TRIAGEoutput.condensed <- TRIAGEoutput[TRIAGEoutput$TRIAGEhit == "Yes", c("EntrezID", "GeneSymbol", "ConfidenceCategory", "TRIAGEhit", "Pathway", "InteractingGenes", "NetworkGenePathways")]

        ########################
        ######## Pathway Output
        #########################
        FinalEnrichment.df <- pathEnrich
        
        
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
        downloadDir <- paste0(outDir, "/", "TRIAGEfilesToDownload")
        dir.create(downloadDir)
        setwd(downloadDir)
        
        TRIAGE.cond.output.name <- paste0(inputFilePrefix, "_", "TRIAGEhits.csv")
        Enrichment.cond.output.name <- paste0(inputFilePrefix, "_", "TRIAGEenrichment.csv")
        
        write.csv(TRIAGEoutput.condensed, file = TRIAGE.cond.output.name)
        write.csv(FinalEnrichment.condensed, file = Enrichment.cond.output.name)

      ######################
      ## Switch to 'Enriched Pathways' tab and display partial results
      observe({
        if(completed) {
          updateTabsetPanel(session, "inTabset", selected = "enrichedPathways")

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
              sprintf('<script>function extLink() {
                          alert("You are leaving the NIH website! This external link provides additional information that is consistent with the intended purpose of this site. NIH cannot attest to the accuracy of a non-federal site. Linking to a non-federal site does not constitute an endoresment by NIH or any of its employees of the sponsors or the information and products presented on the site. You will be subject to the destination site privacy policy when you follow the link.");
                       }</script>                        
                       <form target="_blank" enctype="multipart/form-data" method="post" action="http://www.genome.jp/kegg-bin/mcolor_pathway">
                       <input type="hidden" name="map" value="%s0%s">
                       <input type="hidden" name="unclassified" value="%s">
                       <input type="hidden" name="s_sample" value="color">
                       <input type="hidden" name="mode" value="color">
                       <input type="hidden" name="reference" value="white">
                       <input type="submit" onclick="extLink()" style="font-face: \'Comic Sans MS\'; font-size: larger; color: teal; background-color: powderblue; border: 0 none;"value="%s"></form>', organismAbbr, pathwayID, myGeneLabels, pathwayName)
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
              myBlueGene <- "" 
              myBlueGeneID <- ""
              myBlueGeneLabel <- ""
              myBlueGeneLabels <- ""
              myBlueGeneIDs <- ""
              myRedGene <- ""
              myRedGeneID <- ""
              myRedGeneLabel <- ""
              myRedGeneLabels <- ""
              myRedGeneIDs <- ""


              for(j in 1:length(myGenes))
              {
                # Color the genes on the original hit list BLUE
                if(grepl(myGenes[j], myOriginalHits) ){
                  myBlueGene <- paste(myBlueGene, fontBlue(myGenes[j]), sep = ",")
                  myBlueGeneID <- as.character(siRNA.Score$EntrezID[which(siRNA.Score$GeneSymbol == myGenes[j])])
                  myBlueGeneLabel <- capture.output(cat(myBlueGeneID, "\t#abebc6,blue\t#abebc6,blue"))
                  myBlueGeneLabels <- paste(myBlueGeneLabels, myBlueGeneLabel, sep="\n")
                  myBlueGeneIDs <- capture.output(cat(myBlueGeneIDs, myBlueGeneLabel))
                }else{
                  # Color the genes on the original hit list RED
                  myRedGene <- paste(myRedGene, fontRed(myGenes[j]), sep = ",")
                  myRedGeneID <- as.character(siRNA.Score$EntrezID[which(siRNA.Score$GeneSymbol == myGenes[j])])
                  myRedGeneLabel <- capture.output(cat(myRedGeneID, "\t#ddccff,red\t#ddccff,red"))
                  myRedGeneLabels <- paste(myRedGeneLabels, myRedGeneLabel, "\n")
                  myRedGeneIDs <- capture.output(cat(myRedGeneIDs, myRedGeneLabel))
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
            # Chang column name from 'Genes' to 'TotalGenes'
            colnames(pathEnrich)[which(names(pathEnrich) == "Genes")] <- "TotalGenes"
            return(pathEnrich)
            completed2 <- TRUE
          }, escape = FALSE)
        }
      })

      # Create the 'Gene List' tab
      # output$geneHits <- renderUI({
      #   updateTabsetPanel(session, "inTabset", selected = "geneHits")
        
        output$triageHits <- renderDataTable({
          dat <- datatable(TRIAGEoutput.condensed, rownames = FALSE, options = list(paging=T, autoWidth = F, scrollX = F
                                                                                    , columnDefs = list(list(width = '200px'
                                                                                                             , length = '400px'
                                                                                                             , targets = c(4,5,6)
                                                                                                             ,render = JS(
                                                                                                               "function(data, type, row, meta) {"
                                                                                                               ,"return type === 'display' && typeof data === 'string' && data.length > 30 ?"
                                                                                                               ,"'<span title=\"' + data + '\">' + data.substr(0, 25) + '...</span>' : data;"
                                                                                                               ,"}"))))) 
          return(dat)
        })
        
        output$geneList <- renderDataTable({
            EnrichColumns.index <- NULL
            TRIAGErenamed <- TRIAGEoutput
            
             
            for (i in 1:iterationNum){
              pathColumn <- which(colnames(TRIAGErenamed) == paste0("KEGG.class.iteration", i))
              netColumn <- which(colnames(TRIAGErenamed) == paste0("Network.class.iteration", i))
              colnames(TRIAGErenamed)[pathColumn] <- paste0("PathwayEnrichment", i)
              colnames(TRIAGErenamed)[netColumn] <- paste0("NetworkEnrichment", i)
              EnrichColumns.index <- c(EnrichColumns.index, pathColumn, netColumn)
            }
            
            TRIAGEiterations <- TRIAGErenamed[, c(1:(min(EnrichColumns.index)-2), EnrichColumns.index, (max(EnrichColumns.index)+1):length(TRIAGEoutput))]
            Iteration.index <- c(which(colnames(TRIAGEiterations) == cutoffType), grep('NetworkEnrichment', names(TRIAGEiterations)))
            
            
            totalRow <- data.frame(matrix(NA,1,length(TRIAGEiterations)))
            colnames(totalRow) <- colnames(TRIAGEiterations)
            hitsDataFrame <<- data.frame(matrix(0, length(Iteration.index), 4))
            colnames(hitsDataFrame) <- c('Iteration', 'Total', 'High-conf', 'Med-conf')
            
            for (l in 1:length(Iteration.index)){
              totalHits <- length(which(TRIAGEiterations[Iteration.index[l]] == cutoffHigh))
              totalHighConf <- length(which(TRIAGEiterations[Iteration.index[l]] == cutoffHigh & TRIAGEiterations$ConfidenceCategory == "HighConf"))
              totalMedConf <- length(which(TRIAGEiterations[Iteration.index[l]] == cutoffHigh & TRIAGEiterations$ConfidenceCategory == "MedConf"))
              totalRow[1, Iteration.index[l]] <- totalHits
              hitsDataFrame[l, ] <- c(l - 1, totalHits, totalHighConf, totalMedConf)
            }
            
            totalRow[1,1] <- "Total"
            TRIAGEiterations <- rbind(totalRow, TRIAGEiterations)
            
            # View the dataframe 
            geneHitsToPlot <<- data.frame(hitsDataFrame)
            
            # Highlight the 'Total' row using formatStyle()
            dat <- datatable(TRIAGEiterations, rownames = FALSE, options = list(paging=T, autoWidth = F, scrollX = F, 
                                                                                columnDefs = list(list(width = '200px', length = '400px',
                                                                                                         targets = c((length(TRIAGEiterations)-2), (length(TRIAGEiterations)-1)),
                                                                                                         render = JS(
                                                                                                           "function(data, type, row, meta) {",
                                                                                                           "return type === 'display' && typeof data === 'string' && data.length > 30 ?",
                                                                                                           "'<span title=\"' + data + '\">' + data.substr(0, 25) + '...</span>' : data;",
                                                                                                           "}"))))) %>%
              formatStyle('GeneSymbol', target = 'row', backgroundColor = styleEqual(c('Total'), c('orange')))
            return(dat)
          })
          #message("completed geneList tab")
               
          # Create plots showing the numbers of gene hits by iteration
          output$geneHitsByIteration <- renderPlot({
            geneHitsToPlot.melt <- melt(geneHitsToPlot, id.vars = "Iteration")
            
            ggplot(data = geneHitsToPlot.melt, aes(x = as.numeric(Iteration), y = as.numeric(value), group = variable, color = variable)) +
              geom_line() + geom_point() + labs(x = "Enrichment Iteration", y = "Number of Gene Hits") + theme_light() +
              scale_colour_discrete("") + scale_shape_manual("") + 
              annotation_custom(tableGrob(geneHitsToPlot, rows=NULL), xmin=2, xmax=iterationNum, ymin=numTotal/2, ymax=numTotal) + 
              theme(
                axis.text=element_text(size=12),
                axis.title=element_text(size=14,face="bold")
              )
            
            #theme(axis.title.x = element_text(size = rel(1.8)) + theme(axis.title.y = element_text(size = rel(1.8))
            #ggtitle("Gene Enrichment By Iteration") + theme_bw() + theme(plot.title = element_text(hjust=0.5))
          })
          
          output$pathwayEnrich.cond <- renderDataTable({
            dat <- datatable(FinalEnrichment.condensed, rownames = FALSE, options = list(paging=T, autoWidth = F, scrollX = F
                                                                                      , columnDefs = list(list(width = '200px'
                                                                                                               , length = '400px'
                                                                                                               , targets = c(7,8)
                                                                                                               ,render = JS(
                                                                                                                 "function(data, type, row, meta) {"
                                                                                                                 ,"return type === 'display' && typeof data === 'string' && data.length > 30 ?"
                                                                                                                 ,"'<span title=\"' + data + '\">' + data.substr(0, 25) + '...</span>' : data;"
                                                                                                                 ,"}"))))) #%>%
              #formatStyle("EnrichScore", background = styleColorBar(c(0,1), 'lightblue')) 
            return(dat)
          })
          #})
      
      
      
      ################################
      # Create the 'Network Graph' tab
      ################################    
      output$myNetworkGraph <- renderUI({
        message("Inside NetworkGraph")
        updateTabsetPanel(session, "inTabset", selected = "myNetworkGraph")

        colnames(sigPathways) <- c("Pathways", "TotalGenes", "HitGenes")

        shinyInput <- function(FUN,id,num,...) {
          inputs <- character(num)
          for (i in seq_len(num)) {
            inputs[i] <- as.character(FUN(paste0(id,i),label=NULL,...))
          }
          inputs
        }

        # Select rows (pathways) to generate network graph
        rowSelect <- reactive({

          rows=names(input)[grepl(pattern = "srows_",names(input))]
          paste(unlist(lapply(rows,function(i){
            if(input[[i]]==T){
              return(substr(i,gregexpr(pattern = "_",i)[[1]]+1,nchar(i)))
            }
          })))
        })

        # Show selected row (pathways)
        observe({
          updateTextInput(session, "mySelection", value = rowSelect() ,label = "Your selected pathway IDs:" )
          selectedRows <<- rowSelect()
        })

        output$myNetworkGraph = renderDataTable({
          datatable(cbind(Pick=shinyInput(checkboxInput,"srows_",nrow(sigPathways),value=NULL,width=1), sigPathways[, colnames(sigPathways), drop=FALSE]),
                  options = list(orderClasses = TRUE,
                                 ordering = F,
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

        ## Generating network graph
        observeEvent(input$submitButton, {
          # Only select 1-3 pathways
          if(length(selectedRows) == 0){
            showNotification("You must select at least ONE pathway!", duration = 5,  type=c("error"))
          }
          else if(length(selectedRows) > 3){
            showNotification("You selected more than 3 pathways! Please remove the extra pathways before your re-submission", duration = 10, type=c("error"))
          }
          else{
            ## Show progress bar
            # Create a Progress object
            progress1 <- shiny::Progress$new()
            # Make sure it closes when we exit this reactive, even if there's an error
            on.exit(progress1$close())
            progress1$set(message = "Generating network graph....", value = 0.3)
            
            message(selectedRows)
            source(paste0(scriptDir, "Ranking_plusComments_v3.R"), local = TRUE)
            progress1$inc(1/2)
            Generate_NetworkGraph(selectedRows, organism)
            progress1$inc(1)
            
            #############
            ## PathNet ##
            #############
            saveEdgebundle(Chimera1, paste0(PathNetName.output, "1Degree.html"), selfcontained = TRUE)
            saveEdgebundle(Chimera2, paste0(PathNetName.output, "2Degree.html"), selfcontained = TRUE)
            
            # Remove existing html file in www folder
            Sys.chmod(wwwDir, mode = "0777")
            
            ## Display 1st dimension
            # legend
            output$graphLegend1 <- renderUI({
              HTML(graphLegend)
            })
            #output$graphView1i <- renderEdgebundle({
            output$graphView1i <- renderUI({
              
              file.copy(paste0(PathNetName.output, "1Degree.html"), paste0(wwwDir, paste0(PathNetName.output, "1Degree.html")), overwrite = TRUE)
              tags$iframe(
                seamless="seamless",
                src=paste0(PathNetName.output, "1Degree.html"),
                scrolling = 'no',
                height=700, 
                width=700
              )
            })
            
            # Display 2nd dimension)
            # legend
            output$graphLegend2 <- renderUI({
              HTML(graphLegend)
            })
            
            output$graphView2i <- renderUI({
              # Copy HTML files to www directory for display in iframe
              file.copy(paste0(PathNetName.output, "2Degree.html"), paste0(wwwDir, paste0(PathNetName.output, "2Degree.html")), overwrite = TRUE)
              tags$iframe(
                seamless="seamless",
                src=paste0(PathNetName.output, "2Degree.html"),
                scrolling = 'no',
                height=700, 
                width=700
              )              
            })
            
            output$PathNetTable <- renderDataTable({
              dat <- datatable(Scores_nodes_and_edges, rownames = FALSE, options = list(paging=T, autoWidth = F, scrollX = F
                                                                                           , columnDefs = list(list(width = '200px'
                                                                                                                    , length = '400px'
                                                                                                                    , targets = c(4:(length(Scores_nodes_and_edges)-1))
                                                                                                                    ,render = JS(
                                                                                                                      "function(data, type, row, meta) {"
                                                                                                                      ,"return type === 'display' && typeof data === 'string' && data.length > 30 ?"
                                                                                                                      ,"'<span title=\"' + data + '\">' + data.substr(0, 25) + '...</span>' : data;"
                                                                                                                      ,"}")))))
              return(dat)
            })
            
            ################
            ## D3 Network ##
            ################
            # Display in D3 (1st dimension)
            output$networkView1D3 <- renderForceNetwork({
              # build the two data frames - 'links' and 'nodes'
              g11_links <- g11_d3$links
              g11.links.original <<- g11_links
              g11_nodes <- NodeInfo1

              # delete unneeded data columns
              g11_nodes$ID <- NULL
              g11_nodes$GeneMappingID <- NULL
              # add column names
              colnames(g11_nodes) <- c("name", "group")
              g11.nodes.original <<- g11_nodes

              # Add pathway name as new nodes (with spacd removed - gsub) 
              # and get the indices of the pathway nodes to update the corresponding links
              if(length(selectedRows) == 3){
                g11_nodes[nrow(g11_nodes) + 1,] = c(paste0('[Pathway]', path1_name), 1)
                index1 <<- nrow(g11_nodes)
                g11_nodes[nrow(g11_nodes) + 1,] = c(paste0('[Pathway]', path2_name), 3)
                index2 <<- nrow(g11_nodes)
                g11_nodes[nrow(g11_nodes) + 1,] = c(paste0('[Pathway]', path3_name), 2)
                index3 <<- nrow(g11_nodes)
                # Add Nodesize column
                g11_nodes$nodesize <- 1
                g11_nodes[index1, 'nodesize'] <- 20
                g11_nodes[index2, 'nodesize'] <- 20
                g11_nodes[index3, 'nodesize'] <- 20
              }else if(length(selectedRows) == 2){
                g11_nodes[nrow(g11_nodes) + 1,] = c(paste0('[Pathway]', path1_name), 1)
                index1 <<- nrow(g11_nodes)
                g11_nodes[nrow(g11_nodes) + 1,] = c(paste0('[Pathway]', path2_name), 3)
                index2 <<- nrow(g11_nodes)
                # Add Nodesize column
                g11_nodes$nodesize <- 1
                g11_nodes[index1, 'nodesize'] <- 20
                g11_nodes[index2, 'nodesize'] <- 20
              }else if(length(selectedRows) == 1){
                g11_nodes[nrow(g11_nodes) + 1,] = c(paste0('[Pathway]', path1_name), 1)
                index1 <<- nrow(g11_nodes)
                # Add Nodesize column
                g11_nodes$nodesize <- 1
                g11_nodes[index1, 'nodesize'] <- 20
              }

              # # Update the g11_links with the added pathway nodes
              for (i in 1:(index1 - 1)){
                if(g11_nodes[i, 'group'] == 1){
                  g11_links[nrow(g11_links) + 1,] <- c((index1-1), (i-1))
                }else if(g11_nodes[i, 'group'] == 3){
                  g11_links[nrow(g11_links) + 1,] <- c((index3-1), (i-1))
                }else if(g11_nodes[i, 'group'] == 2){
                  g11_links[nrow(g11_links) + 1,] <- c((index2-1), (i-1))
                }
              }
              
              # Create a copy of network nodes and links for visNetwork
              visNodes11 <<- g11_nodes
              visLinks11 <<- g11_links

              # Save to a html file
              forceNetwork(Links = g11_links, Nodes = g11_nodes,
                           Source = 'source', Target = 'target', NodeID = 'name', Nodesize = 'nodesize',
                           Group = 'group', opacity = 0.6, bounded = FALSE, opacityNoHover = TRUE, fontSize = 10, zoom = TRUE) %>%
                saveNetwork(file = paste0(PathNetName.output, "1Degree.D3.html"))
              
              # Create force directed network plot
              forceNetwork(Links = g11_links, Nodes = g11_nodes,
                           Source = 'source', Target = 'target', NodeID = 'name', Nodesize = 'nodesize',
                           Group = 'group', opacity = 0.6, bounded = FALSE, opacityNoHover = TRUE, fontSize = 10, zoom = TRUE)

              # Add a search box for D3 network
            })
                        
            ## Display in D3 (2nd dimension)
            output$networkView2D3 <- renderForceNetwork({
              # build the two data frames - 'links' and 'nodes'
              g22_links <- g22_d3$links
              g22.links.original <<- g22_links
              g22_nodes <- NodeInfo2
              # delete unneeded data columns
              g22_nodes$ID <- NULL
              g22_nodes$GeneMappingID <- NULL
              # add column names
              colnames(g22_nodes) <- c("name", "group")
              g22.nodes.original <<- g22_nodes
              
              # Add pathway nodes and get the indices for the pathway nodes
              if(length(selectedRows) == 3){
                g22_nodes[nrow(g22_nodes) + 1,] = c(paste0('[Pathway] ', path1_name), 1)
                index1 <<- nrow(g22_nodes)
                g22_nodes[nrow(g22_nodes) + 1,] = c(paste0('[Pathway] ', path2_name), 3)
                index2 <<- nrow(g22_nodes)
                g22_nodes[nrow(g22_nodes) + 1,] = c(paste0('[Pathway] ', path3_name), 2)
                index3 <<- nrow(g22_nodes)
                # Add Nodesize column
                g22_nodes$nodesize <- 1
                g22_nodes[index1, 'nodesize'] <- 20
                g22_nodes[index2, 'nodesize'] <- 20
                g22_nodes[index3, 'nodesize'] <- 20
              }else if(length(selectedRows) == 2){
                g22_nodes[nrow(g22_nodes) + 1,] = c(paste0('[Pathway] ', path1_name), 1)
                index1 <<- nrow(g22_nodes)
                g22_nodes[nrow(g22_nodes) + 1,] = c(paste0('[Pathway] ', path2_name), 3)
                index2 <<- nrow(g22_nodes)
                # Add Nodesize column
                g22_nodes$nodesize <- 1
                g22_nodes[index1, 'nodesize'] <- 20
                g22_nodes[index2, 'nodesize'] <- 20
              }else if(length(selectedRows) == 1){
                g22_nodes[nrow(g22_nodes) + 1,] = c(paste0('[Pathway] ', path1_name), 1)
                index1 <<- nrow(g22_nodes)
                # Add Nodesize column
                g22_nodes$nodesize <- 1
                g22_nodes[index1, 'nodesize'] <- 20
              }
              
              # Update the g11_links with the added pathway nodes
              for (i in 1:(index1 - 1)){
                if(g22_nodes[i, 'group'] == 1){
                  g22_links[nrow(g22_links) + 1,] <- c((index1-1), (i-1))
                }else if(g22_nodes[i, 'group'] == 3){
                  g22_links[nrow(g22_links) + 1,] <- c((index3-1), (i-1))
                }else if(g22_nodes[i, 'group'] == 2){
                  g22_links[nrow(g22_links) + 1,] <- c((index2-1), (i-1))
                }
              }      
              
              # Create a copy of network nodes and links for visNetwork
              visNodes22 <<- g22_nodes
              visLinks22 <<- g22_links

              # Save to a html file
              forceNetwork(Links = g22_links, Nodes = g22_nodes,
                           Source = 'source', Target = 'target', NodeID = 'name', Nodesize = 'nodesize',
                           Group = 'group', opacity = 0.6, bounded = FALSE, opacityNoHover = TRUE, fontSize =10, zoom = TRUE) %>%                            
                saveNetwork(file = paste0(PathNetName.output, "2Degree.D3.html"))
              
              # Create force directed network plot
              forceNetwork(Links = g22_links, Nodes = g22_nodes,
                           Source = 'source', Target = 'target', NodeID = 'name', Nodesize = 'nodesize',
                           Group = 'group', opacity = 0.6, bounded = FALSE, opacityNoHover = TRUE, fontSize =10, zoom = TRUE)
              
              # Add a search box for D3 network
            })
            
            
            ################
            ## visNetwork ##
            ################
            
            # output$networkView3vis <- renderVisNetwork({
            #   # Create force directed network plot
            #   visNetwork(g11_vis$nodes, g11_vis$edges)
            # })
            # output$networkView4vis <- renderVisNetwork({
            #   # Create force directed network plot
            #   visNetwork(g22_vis$nodes, g22_vis$edges)
            # })
          
            # Switch tab
            updateTabsetPanel(session, "inTabset", selected = "networkViews")
            ########################################        

            # Switch tab
            updateTabsetPanel(session, "inTabset", selected = "graphViews")
          }
        })
      })
      message("Outside NetworkGraph")
    
      ## Create the 'Download' tab
      output$downloadFiles <- renderUI({
        updateTabsetPanel(session, "inTabset", selected = "downloads")

        outputFiles <- list.files(path = './')
        out <- c("<br><b>All files from your TRIAGE analysis for download:</b><br>")

        for(i in 1:length(outputFiles)){
          out <- paste(out, outputFiles[i], sep = "<br>")
        }
        out <- paste0(out, "<br><br>")
        HTML(out)
      })

      ## Download all output data
      output$downloadButton <- downloadHandler(

        filename = function(){
          paste("TRIAGE_analysis_output", "zip", sep=".")
        },
        content = function(filename){
          outputFiles <- list.files(path = './')
          zip(zipfile=filename, files = outputFiles)
        },
        contentType = "application/zip"
      )

      message("Download content completed")
    })

    # Contact us by sending email to triage@niaid.nih.gov
    # Send email using mailR package
    # in order for this to work, set the java path for R on commandline in a terminal
    # $sudo ln -f -s $(/usr/libexec/java_home)/jre/lib/server/libjvm.dylib /usr/local/lib

    output$contactUS <- renderUI({

      tagList(
            textInput("userName", "Your Name:", placeholder = "your name"),
            textInput("from", "Your Email Address:", placeholder = "your email address"),
            textInput("subject", "Subject:", placeholder = "Subject"),
            textAreaInput(inputId = "message", label= "Your Email Content:", width = "600px", height = "200px", resize = "vertical", placeholder = "Enter your message here"),
            checkboxInput("contact_not_a_robot", "I'm not a robot*", value = FALSE),
            actionButton("send", " Send email", icon("send-o"), style="padding:4px; font-size:120%; color: #fff; background-color: rgb(1, 81, 154); border-color: #2e6da4")
      )
    })

    # Send
    observeEvent(input$send,{

      # Send email if the 'Send email!' button is clicked and the 'I am not a robot' checked
      if( is.null(input$send) || input$send==0 || !input$contact_not_a_robot){
        return(NULL)
      }

      isolate({
        # Send the email to TRIAGE team
        send.mail(from = input$from,
                  to <- c("jian.song@nih.gov","sakatz@nih.gov"),
                  #to <- c("jian.song@nih.gov"),
                  subject = input$subject,
                  body = paste(paste0("This email is from: ", input$userName, " [", input$from, "]"), input$message, sep="<br/>"),
                  smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = "triage.lisb@gmail.com", passwd = "LISB@NIH", ssl = TRUE),
                  authenticate = TRUE,
                  html = TRUE,
                  send = TRUE)

        # Send email to the user
        send.mail(from = input$from,
                  to = input$from,
                  subject = input$subject,
                  body = paste("Your email was sent!", input$message, sep="<br/>"),
                  smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = "triage.lisb@gmail.com", passwd = "LISB@NIH", ssl = TRUE),
                  authenticate = TRUE,
                  html = TRUE,
                  send = TRUE)
      })

      # Replace the contact me form page with a message after sending the email
      output$contactUS <- renderText({
        "Your email was sent!"
      })
    })

    ## User documentation
    output$documentation <- renderUI({
      tags$iframe(
        seamless="seamless",
        src="UserGuide_V1.pdf",
        height=700,
        width=800
      )
    })
    
    ## Change log
    output$changeLog <- renderUI({
      HTML("Change Log - version and update history")
    })
    
    # Set this to "force" instead of TRUE for testing locally (without Shiny Server)
    session$allowReconnect(TRUE)
  }

  #####################
  # Run the application
  shinyApp(ui = ui, server = server)

#}
