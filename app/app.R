# TRIAGE app
# edits March 14, 2019
#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
options(scipen = 999)
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
library(stringr)
Sys.setenv(R_ZIPCMD="/usr/bin/zip")

# used to set home directory for local development. 
# Assign the home_string for folder where TRIAGE has been downloaded to
#
# home_string = '/Users/username/Documents/'
# Sys.setenv(HOME = home_string)
# setwd('~')

#setting
#override scientific notation to avoid numeric mis assignments
options(scipen = 999)

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

    ##################################################
    # Define UI for application
    ui <- fluidPage(
      
      # Capture user access information
      tags$head(
        tags$title("TRIAGE - Throughput Ranking by Iterative Analysis of Genomic Enrichment"),
        tags$script(src="getIP.js"),
        # sources below are for d3 network graph layout and interaction
        tags$script(src="http://mbostock.github.io/d3/talk/20111116/d3/d3.js"),
        tags$script(src="http://mbostock.github.io/d3/talk/20111116/d3/d3.layout.js"),
        tags$script(src="custom_network.js"),
        tags$script(src="custom_network2.js")
      ),

      # style
      theme = "./css/triage.css",

      # use shinyjs
      useShinyjs(), br(),

      # Sidebar with a slider input for number of bins
      sidebarLayout(

        sidebarPanel(
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
          # textInput("cutoff_valueH", "High-conf Cutoff Value", placeholder = "High-conf cutoff"),
          textInput("cutoff_valueH", "High-conf Cutoff Value", value = "1"),
          bsPopover("cutoff_valueH", "High confidence cutoff value:", "Please enter a value for high confience cutoff, use \"-\" sign for negative value", placement = "bottom", trigger = "hover", options = NULL),
          # textInput("cutoff_valueM", "Med-conf Cutoff Value", placeholder = "Med-conf cutoff"),
          textInput("cutoff_valueM", "Med-conf Cutoff Value", value = "0.5"),
          bsPopover("cutoff_valueM", "Medium confidence cutoff value:", "Please enter a value for medium confience cutoff, use \"-\" sign for negative value", placement = "bottom", trigger = "hover", options = NULL),
          checkboxInput("includeBackground", "Add genome background"),
          bsPopover("includeBackground", "To include known coding genes that are not on your input gene list as background", placement = "bottom", trigger = "hover", options = NULL),          actionButton("goButton", "Analyze my data",
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
                                          dataTableOutput("geneHitsTableByIteration"),
                                          plotOutput("geneHitsByIteration")),
                                 tabPanel(title = "Pathway Enrichments", value = "pathwayEnrich.cond",
                                          dataTableOutput("pathwayEnrich.cond"))
                     )
            ),
            tabPanel(title = "Network", value = "myNetworkGraph",
                     h4('Please select your (1-3) pathways for network graph analysis'), hr(),
                     div(style="display:inline-block",textInput(inputId="mySelection", label="Your selected pathway IDs", value = 0.0)),
                     div(style="display:inline-block",uiOutput("submitGraph")),
                     div(style="display:inline-block",uiOutput("link2Graph")),
                     dataTableOutput("myNetworkGraph")
            ),
            tabPanel(title = "PathNet", value = "graphViews",
                htmlOutput("spacer4"),
                tabsetPanel(id = 'igraphViews',
                      # Display in igraph
                      tabPanel(title="1st Degree Network", value="graphView1",
                               HTML("<div id='graphView1'></div>"),
                               #htmlOutput("graphLegend1"),
                               htmlOutput("graphView1i", width = "100%", height = "700px")
                      ),
                      tabPanel(title="2nd Degree Network", value="graphView2",
                               HTML("<div id='graphView2'></div>"),
                               #htmlOutput("graphLegend2"),
                               htmlOutput("graphView2i", width = "100%", height = "700px")
                      ),
                      tabPanel(title = "PathNet Table", value = "PathNetTable",
                               HTML("<div id='PathNetTable'></div>"),
                               dataTableOutput("PathNetTable")
                      ),
                      tabPanel(title = "Clicked Pathways Table", value = "ClickedDataTable",
                               HTML("<div id='ClickedTable'></div>"),
                               dataTableOutput("ClickedDataTable")
                      )
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
    # Define server logic
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
        dataDir <<- "~/TRIAGE/app/data/"
      }  
      
      # Read in the input fie
      output$contents <- renderDataTable({
        inFile <- input$file1

        if (is.null(inFile))
          return(NULL)

        data <- read.csv(inFile$datapath, stringsAsFactors = FALSE, header=TRUE)

        ## Complete list of 19191 protein-encoding genes in human genome
        ## ftp://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/Mus_musculus.gene_info.gz
        humanGenes <- read.table(file=paste0(dataDir, "HGNC_19191_genes_with_protein_product_EntrezID_geneSymbole_lookup.txt"), sep="\t", header=TRUE)

        ## Complete list of 23504 genes (mRNAs and ncRNAs) in mouse genome
        ## ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_types/gene_with_protein_product.txt
        mouseGenes <- read.table(file=paste0(dataDir, "Mouse_proteinCoding_sourceMouseMine.txt"), sep="\t", header=TRUE)

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

        # Make sure the EntrezIDs are integers
        data$EntrezID <- as.integer(data$EntrezID)
        data <- data[order(data$EntrezID),]

        # Make a copy of the original input data for later use
        siRNA.Score <<- data

        data <- data %>%
          dplyr::select(GeneSymbol, everything())

        # display the input file dimension
        datatable(data, rownames = FALSE, options = list(paging=TRUE))
      })
      
      output$cutoffTypes <- renderUI({
        inFile2 <- input$file1
        if (is.null(inFile2))
          return(NULL)

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

        selectInput("cutoff_type", "Cutoff Type", pulldown_types, selected="Score")
      })

      # reloads the app
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
          siRNA.Score = siRNA.Score[order(siRNA.Score[,'EntrezID'], -siRNA.Score[,input$cutoff_type]),]
          # remove the duplicated EntrezID row that has a smaller value in cutoff_type column 
          siRNA.Score = siRNA.Score[!duplicated(siRNA.Score$EntrezID),]
        }else{
          # sort the dataframe in acending order
          siRNA.Score = siRNA.Score[order(siRNA.Score[,'EntrezID'], siRNA.Score[,input$cutoff_type]),]
          # remove the duplicated EntrezID row that has a larger value in cutoff_type column 
          siRNA.Score = siRNA.Score[!duplicated(siRNA.Score$EntrezID),]
        }

        ## Show progress bar
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        
        progress$set(message = "Performing Enrichment....", value = 0)
        
        startA = Sys.time()
        
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
        organism <<- input$organism
        organismAbbr <- ifelse(grepl("human", tolower(organism)), 'hsa', 'mmu')

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
        
        if(network_ConfidenceCutoff == 700){
          conf.level = ".highConf.igraph"
        }
        else if(network_ConfidenceCutoff == 400){
          conf.level = ".midConf.igraph"
        }
        else{
          conf.level = ".lowConf.igraph"
        }
        
        dir_begin <<- ifelse('SHINY_SERVER_VERSION' %in% env_names, "/srv/shiny-server/data/Networks/String.", "~/TRIAGE/app/data/Networks/String.")
        
        if(input_netowrk == "Experimental & Database"){
          
          # loads igraph for experimental and database source selections
          load(paste0(dir_begin, tolower(organism), ".exp_and_data", conf.level, ".Rdata"))
          G <- get(paste0("String.", tolower(organism), ".exp_and_data", conf.level))
          
        }
        else if(input_netowrk == "Advanced Options"){
          files2load = paste0(dir_begin, tolower(organism), ".", network_InteractionSources, conf.level, ".Rdata")
          vars2load = paste0("String.", tolower(organism), ".", network_InteractionSources, conf.level)
          for (i in 1:length(files2load)) assign(gsub(".*/","",files2load[i]), load(files2load[i]))
          all.G = list()
          for(i in 1:length(vars2load)) all.G[[i]] = get(vars2load[i])
          empty.G.check = which(sapply(sapply(sapply(all.G, edge.attributes), names), is.null))
          count.G.check = which(sapply(all.G, vcount) <= 0)
          all.checks = union(count.G.check, empty.G.check)
          
          datasources = sapply(strsplit(vars2load, '[,.]'), '[[', 3)
          
          if(length(all.checks)>0){
            all.G <- all.G[-union(count.G.check, empty.G.check)]
            datasources = datasources[-union(count.G.check, empty.G.check)]
          }
          
          if(length(all.G)==0){
            showModal(modalDialog(title="Warning:", HTML("<h3><font color=red>Criteria produced empty network. Session will restart.</font><h3>"),
                                  easyClose = TRUE))
            Sys.sleep(5)
            session$reload()
          }
          else if(length(all.G)==1){
            G <- all.G[[1]]
          }
          else{
            
            for(j in 1:length(all.G)){
              names(edge.attributes(all.G[[j]])) = paste0(names(edge.attributes(all.G[[j]])), '_', datasources[j])
            }
            
            temp = all.G[[1]]
            
            for(i in 2:length(all.G)){
              temp = graph.union(temp, all.G[[i]])
            }
            G <<- temp
           
            cols = edge_attr_names(G)[grep("^weights_",edge_attr_names(G))]
            weight_df = list()
            for(i in 1:length(cols)){
              weight_df[[i]] = eval(parse(text=paste0('E(G)$', cols[i])))
              weight_df[[i]][is.na(weight_df[[i]])] = 0
            }
            weight_df = data.frame(weight_df)
            weight_df = apply(weight_df, 2, as.integer)
            E(G)$weights = rowMax(weight_df)
            indx <- max.col(weight_df, ties.method='first')
            cols = edge_attr_names(G)[grep("^datasource_", edge_attr_names(G))]
            E(G)$datasource = datasources[indx]
            message(vertex_attr_names(G))
            message(edge_attr_names(G))
          }
        }
        #Selected_STRINGnetwork.igraph <- G
        message("Networks Loaded")

        #message(networkType)
        use.only.commnected.components <- c('Yes')

        # Remove all files in the outputDir
        unlink(outputDir, recursive = FALSE)
        
        # Create user-specific directory using system time
        userDir <- format(Sys.time(),"%Y%m%d%H%M%S%ms")
        
        outDir <<- paste0(outputDir, userDir)
        dir.create(outDir)
        
        setwd(outDir)
        
        # Set the output file name
        inputFile <- input$file1
        inputFileName <<- inputFile$name
        inputFilePrefix <<- tools::file_path_sans_ext(inputFileName)

        outputFileName <<- paste0(inputFilePrefix, "_", "TRIAGEoutput_ALL.csv")

        # Pathway types
        pathway.type <<- input$pathway
        
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
        
        pathwayData <<- read.csv(file = paste0(dataDir, "Pathways/KEGG2017_", organism, "_", pathway.type, ".csv"))

        proxyScore <- input$cutoff_type
        
        
        #########################################
        ########## Create data tier of gene lists
        #########################################
        cutoffHigh <<- as.numeric(input$cutoff_valueH)
        cutoffMed <<- as.numeric(input$cutoff_valueM)
        cutoffType <<- input$cutoff_type
        #Ensure cutoff column is as numeric
        siRNA.Score <- 
        
        #Depending on the differnce between the high conf cutoff and the mid conf cutoff assign criteria based on "greater than" or "less than"
        if((cutoffHigh - cutoffMed) > 0){
          siRNA.Score <- siRNA.Score %>%
            mutate(ConfidenceCategory = ifelse(get(cutoffType, as.environment(siRNA.Score)) >= cutoffHigh, "HighConf",
                                               ifelse(get(cutoffType, as.environment(siRNA.Score)) >= cutoffMed & get(cutoffType, as.environment(siRNA.Score)) < cutoffHigh,"MedConf",
                                                      "LowConf")))
        }else{
          siRNA.Score <- siRNA.Score %>%
            mutate(ConfidenceCategory = ifelse(get(cutoffType, as.environment(siRNA.Score)) <= cutoffHigh, "HighConf",
                                               ifelse(get(cutoffType, as.environment(siRNA.Score)) <= cutoffMed & get(cutoffType, as.environment(siRNA.Score)) > cutoffHigh,"MedConf",
                                                      "LowConf")))
        }
        
        ## To include background genes 
        includeBackground <- input$includeBackground
        
        if(includeBackground){
          df_background <- data.frame(EntrezID = df_backgroundGenes$EntrezID, GeneSymbol=df_backgroundGenes$GeneSymbol, ConfidenceCategory = rep("Background", nrow(df_backgroundGenes)))
          # combined input data and background data
          myList <- list(siRNA.Score, df_background)
          siRNA.Score <- rbindlist(myList, fill = TRUE)
        }

        cutoffType <<- input$cutoff_type
        proxyScore <- "ConfidenceCategory"
        iteration <- 1
        counter <- TRUE

        # Get a copy of the original list of high-confidence genes
        originalHits <- siRNA.Score$GeneSymbol[siRNA.Score$ConfidenceCategory == "HighConf"]

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
          
          Hits <- siRNA.Score$EntrezID[siRNA.Score[[proxyScore]] == "HighConf"]

          nonHits <- setdiff(siRNA.Score$EntrezID, Hits)

          outPrefix <- paste("KEGG", iteration, sep = "_")

          # 1) Contraction - [Pathway Analysis]
          message(getwd())
          siRNA.Score <- ComputeEnrichment(pathwayData, Hits, nonHits, outPrefix, siRNA.Score, iteration)
          siRNA.Score <- data.frame(siRNA.Score, temp = 0, stringsAsFactors = FALSE)
          kName1 <- paste0("KEGG.class.iteration", iteration)
          kName2 <- paste0("KEGG.", iteration)
          names(siRNA.Score)[names(siRNA.Score) == "temp"] <- kName1
          siRNA.Score[[kName1]][siRNA.Score$KEGG == "Yes" & (siRNA.Score[[proxyScore]] %in% c("MedConf", "HighConf"))] <- "HighConf"
          siRNA.Score[[kName1]][siRNA.Score$KEGG != "Yes" & (siRNA.Score[[proxyScore]] %in% c("MedConf", "HighConf"))] <- "MedConf"
          names(siRNA.Score)[names(siRNA.Score) == "KEGG"] <- kName2

          hit.Genes <- siRNA.Score$EntrezID[siRNA.Score[[kName1]] == "HighConf"]
          myOrignalGenes <- siRNA.Score$GeneSymbol[siRNA.Score[[kName1]] == "HighConf"]

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
          siRNA.Score[[nName2]][siRNA.Score$EntrezID %in% gNames2 & siRNA.Score[[kName1]] %in% c("MedConf", "HighConf")] <- "HighConf"

          
          
          
          #if(iteration != 1 && identical(siRNA.Score[[nName2]], siRNA.Score[[paste0("Network.class.iteration", iteration-1)]])) {
          if((iteration != 1 && identical(siRNA.Score[[nName2]], siRNA.Score[[paste0("Network.class.iteration", iteration-1)]])) 
             || (iteration >= 5 
                 && identical(siRNA.Score[[nName2]], siRNA.Score[[paste0("Network.class.iteration", iteration-2)]]) 
                 && identical(siRNA.Score[[paste0("Network.class.iteration", iteration-1)]], siRNA.Score[[paste0("Network.class.iteration", iteration-3)]])
                 && (length(siRNA.Score$EntrezID[siRNA.Score[[nName2]]== "HighConf"]) > length(siRNA.Score$EntrezID[siRNA.Score[[paste0("Network.class.iteration", iteration-1)]] == "HighConf"])))) { 
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

        samp <- as.data.frame(tempDF[, 2:ncol(tempDF)])
        pathVector <- sapply(seq(nrow(samp)), function(i) unlist(paste(samp[i, which(!is.na(samp[i, ]))], collapse = ", "))) #Switched to comma seperator from " ; " 

        pathDF <- data.frame(GeneSymbol = tempDF$GeneSymbol, Pathway = pathVector, stringsAsFactors = F)

        ## Save the results into output files in the TRIAGEoutputFiles folder
        triage.Out <<- merge(siRNA.Score, pathDF, by = "GeneSymbol", all = T)

        message(getwd(), "#####")

        write.csv(triage.Out, file = outputFileName, row.names = F)
        final_enriched_pathway_file <- paste0(pathway.type, "_TRIAGE_enrichment_final", ".csv")
        write.csv(pathEnrich, file = final_enriched_pathway_file, row.names = F)

        # Pass the pathway list to build networkGraph
        sigPathways <<- pathEnrich[,c("Pathway", "Genes", "HitGenes")]

        completed <- TRUE
     
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
        
        FinalIterationNetworkColumn <- paste0("Network.class.iteration", iterationNum)
        
        TRIAGEoutput <- TRIAGEoutput %>%
          mutate(TRIAGEhit = ifelse(get(FinalIterationNetworkColumn, envir = as.environment(TRIAGEoutput)) == "HighConf", 
                                    "Yes",
                                    ""))
        ################
        # Generate Matrices for TRIAGE Hits, high conf hits, med conf hits
        #Get filtered TRIAGEhits                                           # This is where I put together what is considered a "hit" by TRIAGE (IAM). Any gene that had a score of 1 in the last network analysis step
        
        TRIAGEhits <<- filter(TRIAGEoutput, TRIAGEhit == "Yes")
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
        
        Screen_Genes.for.network.analysis <- intersect(siRNA.Score$EntrezID, V(G)$name)
        
        Graph <- induced.subgraph(G, Screen_Genes.for.network.analysis)
        
        #Screen_Genes.for.network.analysis <- intersect(siRNA.Score$EntrezID,V(G)$name)
        
        #############################################################
        #            Select sub-graphs from hit genes
        #############################################################
        SubGraph <<- Graph
        
        if(length(E(SubGraph))==0 | length(V(SubGraph))==0){
          showModal(modalDialog(title="Warning:", HTML("<h3><font color=red>Criteria produced empty network. Session will restart.</font><h3>"),
                                easyClose = TRUE))
          Sys.sleep(5)
          session$reload()
        }
        
        
        
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
        #GraphEdgesHitNames$weights = round(rowMeans(GraphEdgesHitNames[3:ncol(GraphEdgesHitNames)], na.rm=T), 0)
        GraphEdgesHitNames <- GraphEdgesHitNames[!duplicated(GraphEdgesHitNames), ]
        source <- target <- weights <- rep(NA, nrow(GraphEdgesHitNames))
        tempGenes <- union(GraphEdgesHitNames$from,GraphEdgesHitNames$to)
        for(i in 1:length(tempGenes)){
          source[which(GraphEdgesHitNames$from == tempGenes[i])] <- i-1
          target[which(GraphEdgesHitNames$to == tempGenes[i])] <- i-1
        }
        
        
        
        GraphEdgesHitNumber <- data.table(source,target)
        GraphEdgesHitNumber$weights = GraphEdgesHitNames$weights
        GraphEdgesHitNumber$datasource = GraphEdgesHitNames$datasource
        
        
        
        GraphNodesHit <- data.frame(GeneMappingID = rep(0:(length(tempGenes)-1)), EntrezID = tempGenes)
        
        
        ###############################################################################
        #                 Add back GeneSymbol and Hit Designation to output
        ###############################################################################
        N = ncol(TRIAGEhits)
        # merge(GraphNodesHit, TRIAGEhits[, c("EntrezID", "GeneSymbol", "TRIAGEhit")],
        GraphNodesHit <-  merge(GraphNodesHit, TRIAGEhits[, c(1, 3, 4, N-1, N)],
                               by.x = "EntrezID", by.y = "EntrezID", all.x = T)
        
        #############**************************************************################    # Now getting a data frame for the edges and a data frame for the nodes
        
        ############# EDGE INFO and NODE INFO are sent to GLOBAL ENVIRONMENT to be used by RANKING
        
        EdgeInfo <<- GraphEdgesHitNumber
        NodeInfo <<- GraphNodesHit
        
        #colnames(NodeInfo)[4] = 'Confidence'
        
        
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
        
        #######Add to TRIAGE output file
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
        TRIAGEoutput.condensed <<- TRIAGEoutput[TRIAGEoutput$TRIAGEhit == "Yes", c("EntrezID", "GeneSymbol", "ConfidenceCategory", "TRIAGEhit", "Pathway", "InteractingGenes", "NetworkGenePathways")]

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
        #setwd(downloadDir)
        setwd('TRIAGEfilesToDownload')
        
        TRIAGE.cond.output.name <<- paste0(inputFilePrefix, "_", "TRIAGEhits.csv")
        Enrichment.cond.output.name <- paste0(inputFilePrefix, "_", "TRIAGEenrichment.csv")
        
        write.csv(TRIAGEoutput.condensed, file = TRIAGE.cond.output.name)
        write.csv(FinalEnrichment.condensed, file = Enrichment.cond.output.name)
        # write.csv(triage.Out, file = outputFileName, row.names = F)

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
              sprintf('<a href="https://www.kegg.jp/kegg-bin/show_pathway?%s0%s" target="_blank" class="btn btn-primary">%s</a>', val1, val2, val3)
            }

            # Add a hyperlink to KEGG mapper
            link2KEGGmapper <- function(organismAbbr, pathwayID,  myGeneLabels, pathwayName) {
              # Check the data submitted by the form using http post method
              #sprintf('<form target="_blank" enctype="multipart/form-data" method="post" action="http://localhost/cgi-bin/display_form_data.cgi">
              # Create a form for each datatable row
              sprintf('<script>function extLink() {
                          alert("You are leaving the NIH website! This external link provides additional information that is consistent with the intended purpose of this site. NIH cannot attest to the accuracy of a non-federal site. Linking to a non-federal site does not constitute an endoresment by NIH or any of its employees of the sponsors or the information and products presented on the site. You will be subject to the destination site privacy policy when you follow the link.");
                       }</script>                        
                       <form target="_blank" enctype="multipart/form-data" method="post" action="https://www.kegg.jp/kegg-bin/mcolor_pathway">
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
                if(grepl(paste0(myGenes[j], ","), myOriginalHits) ){
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
            #View(pathEnrich)
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
            Iteration.index <- c(which(colnames(TRIAGEiterations) == "ConfidenceCategory"), grep('NetworkEnrichment', names(TRIAGEiterations)))
            
            
            totalRow <- data.frame(matrix(NA,1,length(TRIAGEiterations)))
            colnames(totalRow) <- colnames(TRIAGEiterations)
            hitsDataFrame <<- data.frame(matrix(0, length(Iteration.index), 4))
            colnames(hitsDataFrame) <- c('Iteration', 'Total', 'High-conf', 'Med-conf')
            
            for (l in 1:length(Iteration.index)){
              totalHits <- length(which(TRIAGEiterations[Iteration.index[l]] == "HighConf"))
              totalHighConf <- length(which(TRIAGEiterations[Iteration.index[l]] == "HighConf" & TRIAGEiterations$ConfidenceCategory == "HighConf"))
              totalMedConf <- length(which(TRIAGEiterations[Iteration.index[l]] == "HighConf" & TRIAGEiterations$ConfidenceCategory == "MedConf"))
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
        
        
          # Show table of gene hits by iterations 
          # It is used to generate the plot below
          output$geneHitsTableByIteration <- renderDataTable({
            dat1 <- datatable(geneHitsToPlot, options = list(paging = FALSE, searching = FALSE, rownames = FALSE))
            return (dat1)
          })
          
          # Create plots showing the numbers of gene hits by iteration
          output$geneHitsByIteration <- renderPlot({
            geneHitsToPlot.melt <- melt(geneHitsToPlot, id.vars = "Iteration")
            
            ggplot(data = geneHitsToPlot.melt, aes(x = as.numeric(Iteration), y = as.numeric(value), group = variable, color = variable)) +
              geom_line() + geom_point() + labs(x = "Enrichment Iteration", y = "Number of Gene Hits") + theme_light() +
              scale_colour_discrete("") + scale_shape_manual("") + 
              #annotation_custom(tableGrob(geneHitsToPlot, rows=NULL), xmin=2, xmax=iterationNum, ymin=numTotal/2, ymax=numTotal) + 
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
              
            return(dat)
          })
      
      
      
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
                                 bFilter = 0,
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
            startB <- Sys.time()
            
            progress1 <- shiny::Progress$new()
            # Make sure it closes when we exit this reactive, even if there's an error
            on.exit(progress1$close())
            progress1$set(message = "Generating network graph....", value = 0.3)
            
            message(selectedRows)
            source(paste0(scriptDir, "config_jsons.R"), local = TRUE)
            # source(paste0(scriptDir, "Ranking_plusComments_v3.R"), local = TRUE)
            source(paste0(scriptDir, "Ranking_source.R"), local = TRUE)
            progress1$inc(1/2)
            
            
            Generate_NetworkGraph(selectedRows, organism, G)
            
            # Writing fully generated network files for download
            write.csv(rbindlist(json_2df), file = paste0(inputFilePrefix, "_", "Second_Degree_Network.csv"))
            write.csv(rbindlist(json_1df), file = paste0(inputFilePrefix, "_", "First_Degree_Network.csv"))
            
            
            print(paste0("Network Generation Time: ", Sys.time() - startB))
            
            
            observe({
              if(input$inTabset == 'downloads'){
                if(length(input$clickedData)>0){
                  clicker <<- data.frame(jsonlite::fromJSON(input$clickedData))
                  clicker = apply(clicker, 2, unlist)
                  #print(clicker)
                  write.csv(clicker, file = paste0(inputFilePrefix, "_", "Clicked_Pathways.csv"))
                }
                else{
                  write.csv(data.frame('Clicked'=NA), file = paste0(inputFilePrefix, "_", "Clicked_Pathways.csv"))
                }
              }
            })
            
            
            observeEvent(input$clickedData, {
              
              clicker <- data.frame(jsonlite::fromJSON(input$clickedData))
              
              output$ClickedDataTable <- renderDataTable({
                if(nrow(clicker)==0){
                  clicker=NULL
                }
                dat <- datatable(data.frame(clicker), rownames = TRUE)
                return(dat)
              })
              
            })
            
            clicked.path.data <- reactive({
              observeEvent(input$clickedData, {
                return(data.frame(jsonlite::fromJSON(input$clickedData)))
              })
            })
            
            
            
            
            
            progress1$inc(1)
            head(E(G))
            
            #############
            ## PathNet ##
            #############
            
            # Remove existing html file in www folder
            Sys.chmod(wwwDir, mode = "0777")
            
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
      
      print(paste0("TRIAGE Analysis Time: ", Sys.time() - startA))
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
        src="UserGuide_V2.pdf",
        height=700,
        width=800
      )
    })
    
    ## Change log
    output$changeLog <- renderUI({
      HTML("Change Log - version and update history")
    })
    
    # This code will be run after the client has disconnected
    # ie, a user session is closed
    # session$onSessionEnded(function() {
    #   # Delete user/session-specific directory after removing all downloadable files
    #   unlink(outDir, recursive = TRUE)
    # })
    
    # Set this to "force" instead of TRUE for testing locally (without Shiny Server)
    session$allowReconnect(TRUE)
  }

  #####################
  # Run the application
  shinyApp(ui = ui, server = server)

