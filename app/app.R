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
library(igraph)
library(edgebundleR)
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

if (interactive()) {
  
    ##################################################
    # Define UI for application that draws a histogram
    ui <- fluidPage( 
      
      # Capture user access information
      tags$head(
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
          # Parameters to be selected
          selectInput(inputId = "organism",
                      label = "Select your organism:",
                      choices = c("Human", "Mouse")
          ),
          selectInput(inputId = "pathway",
                      label = "Select a pathway to enrich:",
                      choices = c("KEGG", "REACTOME") 
          ),
          selectInput(inputId = "network",
                      label = "Select a PPI database to use:",
                      choices = c("STRING high conf", "STRING medium conf") 
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
            tabPanel(title = "Status", value = "status", 
                     tags$pre(id = "status")
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
            tabPanel(title = "Graph View", value = "graphView", 
                     htmlOutput("graphLegend"),
                     edgebundleOutput("graphView", width = "100%", height = "700px")
            ),
            tabPanel(title = "Download", value = "downloads",
                     htmlOutput("downloadFiles"),
                     downloadButton('downloadButton', 'Download all files')
            )
          ), 
          width = 9
        )
      ),
      
      # Get user inputs (user name and email address)
      # bsModal("modalnew", "Info to Access Your Results:", "goButton", size = "small",
      #         textOutput("textnew"),
      #         textInput("userName", "User Name", placeholder = "your name"),
      #         textInput("userEmail", "User Email Address", placeholder = "your email address"),
      #         tags$head(tags$style("#modalnew .modal-footer{ display:none}")),
      #         # check for valid email
      #         conditionalPanel(condition = "input.userName.length > 2
      #                          && input.userEmail.indexOf('@') > -1
      #                          && input.userEmail.indexOf('.') > -1", 
      #                          modalButton("Submit"))
      # ),
      # Show a footer using the header style
      headerPanel(includeHTML("footer.html"))
    )
    
    ##################################################
    # Define server logic required to draw a histogram
    server <- function(session, input, output) {
  
      # Read in the input fie  
      output$contents <- renderDataTable({
        inFile <- input$file1
        
        if (is.null(inFile))
          return(NULL)
        
        data <- read.csv(inFile$datapath)
        
        # display the input file dimension
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

      # These user information (userName and userEmail) will collected after the modal is closed/submitted
      values = reactiveValues(userInfo = "",   ## This is the text that will be displayed
                              modal_closed=F)  ## This prevents the values$userInfo output from updating until the modal is closed/submitted
      
      # Perfrom enrichment
      observeEvent(input$goButton, {   
        
        ## Open the modal when button clicked
        values$modal_closed <- F
        showModal(modalDialog(
          title = "User Info Needed to Access Results",
          size = "s",
          textInput("userName", "User Name", placeholder = "your name"),
          textInput("userEmail", "User Email Address", placeholder = "your email address"),
          
          ## This footer replaces the default "Dismiss" button with 'footer = modalButton("Submit")'
          footer = actionButton("submit_modal",label = "Submit")
        ))
      })
    
      ## This event is triggered by the actionButton inside the modalDialog
      #  It closes the modal, and by setting values$modal_closed <- T, it
      #  triggers values$userInfo to update.
      
      observeEvent(input$submit_modal,{
        values$modal_closed <- T
        
        if(isValidEmail(input$userEmail)){
          removeModal()
        }
        else{
          reset("userEmail")
        }
        
        # Upon job submission, switch to 'status' tab
        message("switching to status tab")
        updateTabsetPanel(session, "inTabset", selected = "status")

        withCallingHandlers({
          shinyjs::html("status", "")
        
        # Global environmental variables
        envs <- Sys.getenv()
        env_names <- names(envs)

        ## Set up scriptDir, inputDir, outputDir depending on whether this is used 
        # as a standalone tool or on AWS webservice
        if('SHINY_APP' %in% env_names){
          scriptDir <- '/srv/shiny-server/Rscripts/'
        }else{
          scriptDir <- "~/TRIAGE/app/Rscripts/"
        }   
        
        if('SHINY_APP' %in% env_names){
          inputDir <- '/srv/shiny-server/inputOutputs/TRIAGEinputFiles/'
        }else{
          inputDir <- "~/TRIAGE/app/inputOutputs/TRIAGEinputFiles/"
        }  
        
        if('SHINY_APP' %in% env_names){
          outputDir <- '/srv/shiny-server/inputOutputs/TRIAGEoutputFiles/'
        }else{
          outputDir <- "~/TRIAGE/app/inputOutputs/TRIAGEoutputFiles/"
        }  
        
        if('SHINY_APP' %in% env_names){
          dataDir <- '/srv/shiny-server/data/'
        }else{
          dataDir <- "~/TRIAGE/app/data/"
        }
        
        # Get organism name from user input
        organism <- input$organism
        organismAbbr <- ifelse(grepl("human", tolower(organism)), 'hsa', 'mmu')
        print(organism)
        
        ## Source other codes depending on whether this is used 
        # as a standalone tool or on AWS webservice 
        if('SHINY_APP' %in% env_names){
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
        # user.log - with user input information
        logDir <- getwd()
        message(logDir, "***")
        userLoginInfo <- c(input$userName, input$userEmail)
        write.table(userLoginInfo, file = paste0(logDir, '/user_login.log'), append = TRUE)
        
        # access.log 
        # Capture user access information
        IP <- reactive({ input$getIP })
        observe({
          cat(capture.output(str(IP()), split=TRUE))
          userAccessInfo <- capture.output(str(IP()), split=TRUE)
          write.table(userAccessInfo, file = paste0(logDir, '/user_access.log'), append = TRUE, quote = TRUE, sep = " ",
                      eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                      col.names = TRUE, qmethod = c("escape", "double"))
        })
        
        # Create uesr-specific directory 
        userDir <- input$userEmail
        userDir <- gsub("\\.", "_", userDir)
        userDir <- gsub("@", "_", userDir)
        setwd(outputDir)
        outDir <- c(getwd(), "/", userDir)
        
        # Remove all files in the outputDir
        unlink(outputDir, recursive = FALSE)
        
        # Create user-specific directory
        if(!dir.exists(userDir)[1]){
          dir.create(userDir)
        }
        else{
          unlink(outDir, recursive = FALSE)
        }
        setwd(userDir)
        
        # Set the output file name
        inputFile <- input$file1
        inputFileName <- inputFile$name
        inputFilePrefix = unlist(strsplit(inputFileName, split='.csv', fixed=TRUE))[1]
        outputFileName <- paste0(inputFilePrefix, "_", networkType, "_TRIGEouput_ALL.csv")
      
        # 1) Seed Pathway Analysis
        # if('SHINY_APP' %in% env_names){
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
        siRNA.Score <- read.csv((input$file1)$datapath, stringsAsFactors = F)
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

        message(getwd(), "#####")
        
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
          else if(length(selectedRows) > 3){
            showNotification("You selected more than 3 pathways! Please remove the extra pathways before your re-submission", duration = 10, type=c("error"))
          }
          else{
            message(selectedRows)
            source(paste0(scriptDir, "Ranking_plusComments_v2.R"), local = TRUE)
            Generate_NetworkGraph(selectedRows, organism)
            
            #observeEvent(input$submitButton, {
              updateTabsetPanel(session, "inTabset", selected = "graphView")
              
              # Display the HTML result in a Shiny tab
              output$graphView <- renderUI({
                message("Inside graphView tab")
                #updateTabsetPanel(session, "inTabset", selected = "graphView")
                
                output$graphLegend <- renderUI({
                  HTML(graphLegend)
                })
                
                output$graphView <- renderEdgebundle({
                  # Display the network graph
                  Chimera
                })
              })
            #})
          }
        })
        
        # observeEvent(input$submitButton, {
        #   if(length(selectedRows) >= 1 & length(selectedRows) <= 3) {
        #     output$link2Graph <- renderUI({
        #       #a(href="file:///Users/songj11/TRIAGE/app/inputOutputs/TRIAGEoutputFiles/Chimera_STRINGHi_MoTNF.hits.html", target="_blank", "Click here to see the resulting network graph")
        #       if(grepl('shiny', outputDir)){
        #         a("Done! Click here to see the resulting network graph",target="_blank",href="/results/Chimera_STRINGHi_MoTNF.hits.html")
        #       }else{
        #         a("Done! Click here to see the resulting network graph",target="_blank",href="http://localhost/Chimera_STRINGHi_MoTNF.hits.html")
        #       }
        #     })
        #   }
        # })
      })
      message("Outside NetworkGraph")
      
      # observeEvent(input$submitButton, {
      #   updateTabsetPanel(session, "inTabset", selected = "graphView")
      #   
      #   # Display the HTML result in a Shiny tab
      #   output$graphView <- renderUI({
      #     message("Inside graphView tab")
      #     #updateTabsetPanel(session, "inTabset", selected = "graphView")
      #     
      #     output$graphLegend <- renderUI({
      #       HTML(graphLegend)
      #     })
      #     
      #     output$graphView <- renderEdgebundle({
      #       # Display the network graph
      #       Chimera
      #     })
      #   })
      # })
      
      # Create the 'Download' tab
      output$downloadFiles <- renderUI({
        updateTabsetPanel(session, "inTabset", selected = "downloads")
        message(outDir, "*****")
        outputFiles <- list.files(path = './')
        out <- c("<br><b>All files from your TRIAGE analysis for download:</b><br>")
        
        for(i in 1:length(outputFiles)){
          out <- paste(out, outputFiles[i], sep = "<br>")
        }
        out <- paste0(out, "<br><br>")
        HTML(out)
      })
      
      # Download all output data
      output$downloadButton <- downloadHandler(
        
        filename = function(){
          paste("TRIAGE_analysis_output","zip",sep=".")
        },
        content = function(filename){
          #outputFiles <- list.files(path = outputDir)
          outputFiles <- list.files(path = outDir)
          zip(zipfile=filename, files = outputFiles)
        },
        contentType = "application/zip"
      )
      
      message("Download content completed")      
    })

    # Set this to "force" instead of TRUE for testing locally (without Shiny Server)
    session$allowReconnect(TRUE)
  }
    
  #####################
  # Run the application 
  shinyApp(ui = ui, server = server)
}