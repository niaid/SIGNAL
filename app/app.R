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
                      label = "Select a pathway to enrich:",
                      choices = c("KEGG")
          ),
          selectInput(inputId = "network",
                      label = "Select a PPI database to use:",
                      choices = c("STRING high conf", "STRING medium conf")
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
          actionButton("goButton", "Analyze my data",
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
            tabPanel(title = "Gene Hits", value = "geneList",
                     htmlOutput("spacer3"),
                     tabsetPanel(id = 'geneList',
                        # Display Gene hits 
                        tabPanel(title="Lists of Gene Hits", value="geneHits",
                         dataTableOutput("geneList")),
                        
                        # Diplay number of gene hits by iteration
                        tabPanel(title="Gene Hits By Iteration", value="geneHitsByIteration",
                                 plotOutput("geneHitsByIteration"))
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
                      tabPanel(title="1st Dimension Network", value="graphView1",
                               htmlOutput("graphLegend1"),
                               htmlOutput("graphView1i", width = "100%", height = "700px")
                      ),
                      tabPanel(title="2nd Dimension Network", value="graphView2",
                               htmlOutput("graphLegend2"),
                               htmlOutput("graphView2i", width = "100%", height = "700px")
                      )
                )
            ),
            tabPanel(title = "NetworkD3", value = "networkViews",
                htmlOutput("spacer5"),
                tabsetPanel(id = 'networkD3Views',
                      ## Display in networkD3
                      tabPanel(title="1st Dimension D3 Network", value="networkView1",
                               forceNetworkOutput("networkView1D3", width = "100%", height = "700px")),

                      tabPanel(title="2nd Dimension D3 Network", value="networkView2",
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
      
      # Read in the input fie
      output$contents <- renderDataTable({
        inFile <- input$file1

        if (is.null(inFile))
          return(NULL)

        data <- read.csv(inFile$datapath, stringsAsFactors = F,header=TRUE)
        
        # # Check for duplicated GeneSymbols
        # if(anyDuplicated(data$GeneSymbol)){
        #   showModal(modalDialog(title="User Input Errors", HTML("<h4><font color=red>Duplicated GeneSymbols were found! <br><br>Please remove the duplicates and reload your input file.</font><h4>")))
        # }        
        # 
        # # Check for duplicated GeneSymbols
        # if(anyDuplicated(data$EntrezID)){
        #   showModal(modalDialog(title="User Input Errors", HTML("<h4><font color=red>Duplicated EntrezIDs were found! <br><br>Please remove the duplicates and reload your input file.</font><h4>")))
        # }        
        
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
          } else if(input$organism == "Mouse"){
            library('org.Mm.eg.db')
            x <- org.Mm.egSYMBOL2EG
          }
          mapped_genes <- mappedkeys(x)
          # Change the GeneSymbols in the input file to UPPERCASE
          #overlappingGenes <- intersect(as.character(as.list(mapped_genes)), as.character(toupper(data$GeneSymbol)))
          overlappingGenes <- intersect(as.character(as.list(mapped_genes)), as.character(data$GeneSymbol))
          xx <- as.list(x[overlappingGenes])
          y <- unlist(xx)
          y <- data.frame(GeneSymbol = names(y), EntrezID = y, row.names = NULL)
          numGeneInInput <- nrow(data)
          # Change the GeneSymbols in the input file to UPPERCASE
          #data$GeneSymbol <- toupper(data$GeneSymbol)
          tempData <- merge(x=data, y=y, by="GeneSymbol")
          data <- tempData
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
          } else if(input$organism == "Mouse"){
            library('org.Mm.eg.db')
            x <- org.Mm.egSYMBOL
          }
          mapped_genes <- mappedkeys(x)
          overlappingGenes <- intersect(mapped_genes,data$EntrezID)
          xx <- as.list(x[!is.na(overlappingGenes)])
          y <- unlist(xx)
          y <- data.frame(GeneSymbol = y, EntrezID = names(y), row.names = NULL)

          tempData <- merge(x=data, y=y, by="EntrezID", all.x = TRUE)
          data <- tempData
          data = data[, c(ncol(data), 1:(ncol(data) - 1))]
          
          rm(tempData,x,y,xx)
          message("Input file has a 'EntrezID' column!")
        }
        
        # Populate GeneSymbolcolumn with EntrezIDs if the corresponding GeneSymbols are not available
        for (i in 1:nrow(data)){
          if(grepl('-', data$GeneSymbol[i])){
            data$GeneSymbol[i] <- data$EntrezID[i]
          }else if(is.na(data$GeneSymbol[i])){
            data$GeneSymbol[i] <- data$EntrezID[i]
          }
        }
        
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
      
        ## Upon job submission, switch to 'status' tab
        # message("switching to status tab")
        # updateTabsetPanel(session, "inTabset", selected = "status")

        #withCallingHandlers({
        #   shinyjs::html("status", "")

        # Global environmental variables
        envs <- Sys.getenv()
        env_names <- names(envs)

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

        if('SHINY_SERVER_VERSION' %in% env_names){
          dataDir <- '/srv/shiny-server/data/'
        }else{
          dataDir <- "~/TRIAGE/app/data/"
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

        # Set the network to be used based on user selection
        network <- input$network

        if(tolower(organism) == 'human'){
          if(grepl("high", network)){
            networkType <- 'hSTRINGppi.hi'
          }
          if(grepl("medium", network)){
            networkType <- 'hSTRINGppi.med'
          }
        }
        else if(tolower(organism) == 'mouse'){
          if(grepl("high", network)){
            networkType <- 'mSTRINGppi.hi'
          }
          if(grepl("medium", network)){
            networkType <- 'mSTRINGppi.med'
          }
        }

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

        outputFileName <- paste0(inputFilePrefix, "_", networkType, "_TRIGEouput_ALL.csv")

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

        pathwayData <- read.csv(file = paste0(dataDir, "Pathways/", pathway.type, organism, ".csv"))

        # Get input file
        # siRNA.Score <- read.csv((input$file1)$datapath, stringsAsFactors = F)
        # Populate GeneSymbolcolumn with EntrezIDs if the corresponding GeneSymbols are not available
        # for (i in 1:nrow(siRNA.Score )){
        #   if(siRNA.Score$GeneSymbol[i] == '-'){
        #     siRNA.Score$GeneSymbol[i] <- siRNA.Score$EntrezID[i]
        #   }
        # }        
        
        proxyScore <- input$cutoff_type
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

          outPrefix <- paste(pathway.type, iteration, sep = "_")

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
        pathVector <- sapply(seq(nrow(samp)), function(i) unlist(paste(samp[i, which(!is.na(samp[i, ]))], collapse = " ; ")))

        pathDF <- data.frame(GeneSymbol = tempDF$GeneSymbol, Pathway = pathVector, stringsAsFactors = F)

        ## Save the results into utput files in the TRIAGEoutputFiles folder
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

      ## Switch to 'Enriched Pathways' tab adn display partial results
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
      output$geneList <- renderUI({
          updateTabsetPanel(session, "inTabset", selected = "geneList")

          output$geneList <- renderDataTable({
            # Select only the hit genes in any iteration
            # number of columns to check based on the number of iterations
            numConditions <- NULL

            for (i in 1:(iteration - 1))
            {
              if (i == (iteration - 1)){
                numConditions <- paste(numConditions, "siRNA.Score$KEGG.class.iteration", i, " > ", as.numeric(input$cutoff_valueM), sep="")
              }
              else{
                numConditions <- paste(numConditions, "siRNA.Score$KEGG.class.iteration", i, " > ", as.numeric(input$cutoff_valueM), " | ", sep="")
              }
            }
            
            #print(numConditions)
            subSet <- subset(siRNA.Score, eval(parse(text = numConditions)))

            # Add a new row of "Total" of hit genes in each iteration
            totalRow <- subSet[1,]

            for (k in 1:length(totalRow)) {

              if(grepl("KEGG.class.iteration", colnames(totalRow)[k])) {

                # Get the subset of the dataframe with column value equal to 1
                subSet1 <- subset(subSet, subSet[[colnames(totalRow)[k]]] == input$cutoff_valueH)
                totalRow[1,colnames(totalRow)[k]] <- length(subSet1[[colnames(totalRow)[k]]])
              }
              else{
                totalRow[1,k] <- NA
              }
            }
            totalRow[1,"GeneSymbol"] <- "Total"

            # The table of gene hits in the enriched pathways by iterations
            newSubset <- rbind(totalRow, subSet)
            # Create a dataframe of geneHits for total, high-conf, med-conf cross multiple iterations of the enrichment
            # The gene hit counts are ONLY from the enriched pathways met with pVal cutoff (0.05)
            #             Total      Hign-Conf    Med-conf
            # Original      x            x           x
            # Iteration1    x            x           x
            # Iteration2    x            x           x
            # ........      x            x           x
            # original input data => siRNA.Score
            # enriched data => newSubset
            # dataframe => geneHitsByIterations 
            cutoffType <- input$cutoff_type
            cutoffHigh <- as.numeric(input$cutoff_valueH)
            cutoffMed <- as.numeric(input$cutoff_valueM)
            numHighConf <- 0
            numMedConf <- 0
            
            # Count the high-conf, med-conf, and total numbers of genes in the input data
            for(i in 1:length(siRNA.Score[,cutoffType])){
              if((!is.na(siRNA.Score[i, cutoffType])) & (siRNA.Score[i, cutoffType] >= cutoffHigh)){
                numHighConf = numHighConf + 1
              }else if((!is.na(siRNA.Score[i, cutoffType])) & (siRNA.Score[i,cutoffType] >= cutoffMed)){
                numMedConf= numMedConf + 1
              }
            }
            numTotal <<- numHighConf + numMedConf
            
            # Create the dataframe
            geneHitsByIterations <<- rbind( c('Total', 'High-conf', 'Med-conf'), c(numTotal, numHighConf, numMedConf))

            # Count the high-conf, med-conf, and total numbers of genes in the enriched pathways
            for(j in 1:iterationNum){
              # Get column names for each iteration
              iterationCol <- paste0("KEGG.class.iteration",j)

              # total number of 'hit' genes in each iteration
              totalCount <- sum(siRNA.Score[,iterationCol] >= cutoffMed, na.rm=TRUE)
              
              # number of 'hit' genes from the 'high-conf' gene set in the input data
              highCount <- sum(siRNA.Score[,iterationCol] >= cutoffHigh, na.rm=TRUE)
              
              # number of 'hit' gene from the 'med-conf' gene set in the input data
              medCount <- sum(siRNA.Score[,iterationCol] < cutoffHigh & siRNA.Score[,iterationCol] >= cutoffMed, na.rm=TRUE)
              
              geneHitsByIterations <- rbind(geneHitsByIterations, c(totalCount, highCount, medCount))
            }
            # Use row#1 as the header 
            colnames(geneHitsByIterations) = geneHitsByIterations[1, ]
            geneHitsByIterations = geneHitsByIterations[-1, ]
            # Use 'Iteration' column with row index number
            #geneHitsByIterations <- cbind(0:(nrow(geneHitsByIterations) - 1), geneHitsByIterations)
            geneHitsByIterations <- cbind(0:(nrow(geneHitsByIterations) - 1), geneHitsByIterations)
            colnames(geneHitsByIterations)[1] <- "Iteration"
            
            # View the dataframe 
            geneHitsToPlot <<- data.frame(geneHitsByIterations)
      
            # Highlight the 'Total' row using formatStyle()
             dat <- datatable(newSubset, rownames = FALSE, options = list(paging=TRUE)) %>%
              formatStyle('GeneSymbol', target = 'row', backgroundColor = styleEqual(c('Total'), c('orange')))
              return(dat)
          })
          message("completed geneList tab")
          
          # Create plots how the numbers of gene hits by iteration
          output$geneHitsByIteration <- renderPlot({
            geneHitsToPlot.melted <- melt(geneHitsToPlot, id.vars = "Iteration")

            ggplot(data = geneHitsToPlot.melted, aes(x = as.numeric(Iteration) - 1, y = as.numeric(value), group = variable, color = variable)) +
              geom_line() + geom_point() + labs(x = "Enrichment Iteration", y = "Number of Gene Hits") + theme_light() +
              scale_colour_discrete("") + scale_shape_manual("") + annotation_custom(tableGrob(geneHitsToPlot, rows=NULL), xmin=2, xmax=iterationNum, ymin=numTotal/3, ymax=numTotal*2/3) + 
              theme(
                axis.text=element_text(size=12),
                axis.title=element_text(size=14,face="bold")
              )
              
              #theme(axis.title.x = element_text(size = rel(1.8)) + theme(axis.title.y = element_text(size = rel(1.8))
              #ggtitle("Gene Enrichment By Iteration") + theme_bw() + theme(plot.title = element_text(hjust=0.5))
          })         
      })

      # Create the 'Network Graph' tab
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
            source(paste0(scriptDir, "Ranking_plusComments_v2.R"), local = TRUE)
            progress1$inc(1/2)
            Generate_NetworkGraph(selectedRows, organism)
            progress1$inc(1)
            
            #############
            ## PathNet ##
            #############
            saveEdgebundle(Chimera1, "Chimera_STRINGHi_against_selectedPathways_1st.hits.html", selfcontained = TRUE)
            saveEdgebundle(Chimera2, "Chimera_STRINGHi_against_selectedPathways_2nd.hits.html", selfcontained = TRUE)
            
            # Remove existing html file in www folder
            Sys.chmod(wwwDir, mode = "0777")
            
            ## Display 1st dimension
            # legend
            output$graphLegend1 <- renderUI({
              HTML(graphLegend)
            })
            #output$graphView1i <- renderEdgebundle({
            output$graphView1i <- renderUI({
              
              file.copy("Chimera_STRINGHi_against_selectedPathways_1st.hits.html", paste0(wwwDir, "Chimera_STRINGHi_against_selectedPathways_1st.hits.html"), overwrite = TRUE)
              tags$iframe(
                seamless="seamless",
                src="Chimera_STRINGHi_against_selectedPathways_1st.hits.html",
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
              file.copy("Chimera_STRINGHi_against_selectedPathways_2nd.hits.html", paste0(wwwDir, "Chimera_STRINGHi_against_selectedPathways_2nd.hits.html"), overwrite = TRUE)
              tags$iframe(
                seamless="seamless",
                src="Chimera_STRINGHi_against_selectedPathways_2nd.hits.html",
                scrolling = 'no',
                height=700, 
                width=700
              )              
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
                  g11_links[nrow(g11_links) + 1,] <- c((index2-1), (i-1))
                }else if(g11_nodes[i, 'group'] == 2){
                  g11_links[nrow(g11_links) + 1,] <- c((index3-1), (i-1))
                }
              }
              
              # Create a copy of network nodes and links for visNetwork
              visNodes11 <<- g11_nodes
              visLinks11 <<- g11_links

              # Save to a html file
              forceNetwork(Links = g11_links, Nodes = g11_nodes,
                           Source = 'source', Target = 'target', NodeID = 'name', Nodesize = 'nodesize',
                           Group = 'group', opacity = 0.6, bounded = FALSE, opacityNoHover = TRUE, fontSize = 10, zoom = TRUE) %>%
                saveNetwork(file = 'Chimera_STRINGHi_against_selectedPathways_1st.D3.html')
              
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
                  g22_links[nrow(g22_links) + 1,] <- c((index2-1), (i-1))
                }else if(g22_nodes[i, 'group'] == 2){
                  g22_links[nrow(g22_links) + 1,] <- c((index3-1), (i-1))
                }
              }      
              
              # Create a copy of network nodes and links for visNetwork
              visNodes22 <<- g22_nodes
              visLinks22 <<- g22_links

              # Save to a html file
              forceNetwork(Links = g22_links, Nodes = g22_nodes,
                           Source = 'source', Target = 'target', NodeID = 'name', Nodesize = 'nodesize',
                           Group = 'group', opacity = 0.6, bounded = FALSE, opacityNoHover = TRUE, fontSize =10, zoom = TRUE) %>%                            
                saveNetwork(file = 'Chimera_STRINGHi_against_selectedPathways_2nd.D3.html')
              
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
