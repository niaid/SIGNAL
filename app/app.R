# SIGNAL app
# edits February, 2021
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
# library(shinyMatrix)
library(readr)
library(dplyr)
library(stringi)
library(DT)
library(data.table)
library(igraph)
library(edgebundleR)
library(shinyAce)
library(rJava)
#library(mailR)
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
library(tools)
library(readxl)
Sys.setenv(R_ZIPCMD="/usr/bin/zip")

### TO RUN LOCALLY ###
# Assign the home_string for folder where SIGNAL has been downloaded to
# home_string = '/Users/username/SIGNAL-folder'
# Sys.setenv(HOME = home_string)
# setwd('~')

### Package check and installation (for local development)
# L = c('shiny', 'shinyjs', 'shinyBS', 'readr', 'dplyr', 'stringi', 'DT', 'data.table',
#       'igraph', 'edgebundleR', 'shinyAce', 'networkD3', 'visNetwork',
#       'AnnotationDbi', 'reshape2', 'ggplot2', 'tidyr', 'gridExtra', 'crosstalk', 
#       'htmltools', 'stringr', 'org.Hs.eg.db', 'org.Mm.eg.db', 'mailR', 'rJava')
# 
# package.check <- lapply(L, FUN = function(x) {
#   if (!require(x, character.only = TRUE)) {
#     if(x != 'org.Hs.eg.db' & x != 'org.Mm.eg.db'){
#       install.packages(x, dependencies = TRUE)
#       library(x, character.only = TRUE)
#     }
#     else{
#       if("BiocManager" %in% installed.packages()){
#         BiocManager::install(x, version = "3.8")
#       }
#       else{
#         install.packages("BiocManager")
#         BiocManager::install(x, version = "3.8")
#       }
#     }
#   }
# })

### Errors with rJava ###
# If running a mac and having trouble loading rJava - follow the steps on this site:
# https://zhiyzuo.github.io/installation-rJava/
# If running a pc and having trouble loading rJava:
# https://www.r-statistics.com/2012/08/how-to-load-the-rjava-package-after-the-error-java_home-cannot-be-determined-from-the-registry/

#setting
#override scientific notation to avoid numeric mis assignments
options(scipen = 999)

# global variables
organism <- NULL
organismAbbr <- NULL
originalHits <- NULL
networkType <- NULL

### SIGNAL veriosn
SIGNAL.version <- "Version 1.0"
SIGNAL.update.date <- "November 20, 2020"

### Database versions
KEGG.version <- "2020"
KEGG.update.date <- "November 20, 2020"

STRING.version <- "11.0"
STRING.update.date <- "November 20, 2020"


# Function to validation email addresses
isValidEmail <- function(x) {
  grepl("\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>", as.character(x),
        ignore.case=TRUE)
}

# Set the maximum input file size to 10Mb
options(shiny.maxRequestSize = 10*1024^2)

#if (interactive()) {
    ##################################################
    # Define UI for application
    ui <- fluidPage(
      
      # Capture user access information
      tags$head(
        tags$title("SIGNAL - Selection by Iterative pathway Group and Network Analysis Looping"),
        tags$script(src="getIP.js"),
        # sources below are for d3 network graph layout and interaction
        tags$script(src="mbostock_d3.js"),  # original file: https://mbostock.github.io/d3/talk/20111116/d3/d3.js
        tags$script(src="mbostock_d3_layout.js"),  # original file: https://mbostock.github.io/d3/talk/20111116/d3/d3.layout.js
        tags$script(src="custom_network.js"),
        tags$script(src="custom_network2.js")
      ),

      # Show a header using the header style
      headerPanel(includeHTML("header.html")),
      
      # style
      theme = "./css/signal.css",

      # use shinyjs
      useShinyjs(), br(),

      # Sidebar with a slider input for number of bins
      sidebarLayout(

        sidebarPanel(
          # Global site tag (gtag.js) - Google Analytics 
          tags$head(
            HTML("<!-- Global site tag (gtag.js) - Google Analytics -->"),
            tags$script(HTML("(function(w,d,s,l,i){w[l]=w[l]||[];w[l].push({'gtm.start':
                            new Date().getTime(),event:'gtm.js'});var f=d.getElementsByTagName(s)[0],
                             j=d.createElement(s),dl=l!='dataLayer'?'&l='+l:'';j.async=true;j.src=
                              'https://www.googletagmanager.com/gtm.js?id='+i+dl;f.parentNode.insertBefore(j,f);
                             })(window,document,'script','dataLayer','GTM-56JHMG');"
                        ))
            
          ),
          tags$body(HTML('<noscript><iframe src="https://www.googletagmanager.com/ns.html?id=GTM-56JHMG"
                            height="0" width="0" style="display:none;visibility:hidden"></iframe></noscript>'
                        )
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
                      choices = c("STRING: Experimental & Database", "Advanced Options")
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
          # textInput("cutoff_valueH", "High Confidence Cutoff Value", placeholder = "High-conf cutoff"),
          textInput("cutoff_valueH", "High Confidence Cutoff Value"),
          bsPopover("cutoff_valueH", "High confidence cutoff value:", "Please enter a value for high confience cutoff, use \"-\" sign for negative value", placement = "bottom", trigger = "hover", options = NULL),
          # textInput("cutoff_valueM", "Med-conf Cutoff Value", placeholder = "Med-conf cutoff"),
          textInput("cutoff_valueM", "Medium Confidence Cutoff Value"),
          bsPopover("cutoff_valueM", "Medium confidence cutoff value:", "Please enter a value for medium confience cutoff, use \"-\" sign for negative value", placement = "bottom", trigger = "hover", options = NULL),
          checkboxInput("secondaryCutoff", "Add an Additional Criteria"),
          uiOutput("secondaryColumn_choice"),
          conditionalPanel(
            condition = "input.secondaryCutoff == 1",
                   fluidRow(
                     column(3,
                            selectInput("secondaryDirection", "Direction:", c('\u2265','\u2264'), selected='\u2265')
                            ),
                     column(6,
                            textInput(inputId="secondaryValue", label="Value")
                            )
                   )
          ),
          bsPopover("secondaryCutoff", "Use the values of another column to set a cutoff that all medium and high confidence hits must meet", placement = "bottom", trigger = "hover", options = NULL),
          checkboxInput("includeBackground", "Add Genome Background"),
          bsPopover("includeBackground", "To include known coding genes that are not on your input gene list as background", placement = "bottom", trigger = "hover", options = NULL),          
          actionButton("goButton", "Analyze my data",
                       style="padding:4px; font-size:120%; color: #fff; background-color: rgb(1, 81, 154); border-color: #2e6da4"),
          actionButton("refresh", "Reset", icon("undo"),
                       style="padding:4px; font-size:120%; color: #fff; background-color: rgb(1, 81, 154); border-color: #2e6da4"),
          width = 3
        ),
        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(id = "inTabset",
                      tabPanel(title = "Home", value = "landingPage",
                               htmlOutput("spacer1"),
                               h3("The SIGNAL platform is designed to facilitate robust hit selection from normalized high-throughput datasets.", align = "center"),
                               HTML('<center><img src="images/SIGNAL_schmematic_Home.png" width="1200"></center>'),
                               br(),
                               h4("For a detailed user guide and sample dataset click on the “User Guide” tab.", align = "center", style = "color:grey"),
                               br(),
                               br(),
                               br(),
                               HTML("<center><p>Our manuscript describing the design, development, and testing of SIGNAL can be found <a href='https://www.cell.com/cell-systems/fulltext/S2405-4712(21)00081-8'>here</a>.</p></center>"),
                               HTML("<center><p>To cite: Katz, S., Song, J., Webb, K. P., Lounsbury, N. W., Bryant, C. E., &amp; Fraser, I. D. C. (2021). SIGNAL: A web-based iterative analysis platform integrating pathway and network approaches optimizes hit selection from genome-scale assays. <em>Cell Systems, 12</em>(4), 338-352.e335. <a href='doi:https://doi.org/10.1016/j.cels.2021.03.001'>doi:https://doi.org/10.1016/j.cels.2021.03.001</a></p></center>")
                               #img(src = "images/SIGNAL_schmematic_Home.png", height="75%", width="75%", align = "center")
                      ),
            tabPanel(title = "Input", value = "contents",
                     htmlOutput("spacer2"),
                     # dataTableOutput("contents"),
                     # br(),
                     # br(),
                     # fluidRow(div(style='display:inline-block; padding-left:0px; padding-right:0px; horizontal-align:top;',
                     #   column(width =  6, 
                     #          #div(style='display:inline-block; padding-left:0px; padding-right:0px; horizontal-align:top;',
                     #          textAreaInput(inputId="manualHits", label="Optional: Enter gene IDs that should be kept \n as high confidence hits throughout the analysis:", width="100%", height="200px", value="", placeholder = "7099, 4615, 51135")
                     #   ),
                     #   column(width = 3,
                     #          div(style=' padding-top:50px;',
                     #          radioButtons("manualIDtypes", h5("Gene ID type:"),
                     #                           choices = list("Entrez ID" = "EntrezID", "Gene Symbol" = "GeneSymbol"), selected = "EntrezID")
                     #          ))))
                     # ,
                     
                     tabsetPanel(id = 'contents',
                                 # Input genes table
                                 tabPanel(title = "Mapped Genes", value = "contents",
                                          dataTableOutput("contents"),
                                          br(),
                                          br(),
                                          fluidRow(div(style='display:inline-block; padding-left:0px; padding-right:0px; horizontal-align:top;',
                                                       column(width =  6,
                                                              #div(style='display:inline-block; padding-left:0px; padding-right:0px; horizontal-align:top;',
                                                              textAreaInput(inputId="manualHits", label="Optional: Enter gene IDs that should be kept \n as high confidence hits throughout the analysis:", width="100%", height="200px", value="", placeholder = "7099, 4615, 51135")
                                                       ),
                                                       column(width = 3,
                                                              div(style=' padding-top:50px;',
                                                                  radioButtons("manualIDtypes", h5("Gene ID type:"),
                                                                               choices = list("Entrez ID" = "EntrezID", "Gene Symbol" = "GeneSymbol"), selected = "EntrezID")
                                                              ))))
                                 ),
                                 tabPanel(title = "Unmapped Genes", value = "contents2",
                                          dataTableOutput("contents2")
                                 )
                     )
                     
                     
                     
                     #textAreaInput(inputId="manualHits", label="Optional: Enter gene IDs that should be kept \n as high confidence hits throughout the analysis:", width="100%", height="200px", value="", placeholder = "7099, 4615, 51135"),
                     # div(style="display:grid",textAreaInput(inputId="manualHits", label="Optional: Enter gene IDs that should be kept \n as high confidence hits throughout the analysis:", width="100%", height="200px", value="", placeholder = "7099, 4615, 51135")),
                     # div(style="display:grid",radioButtons("manualIDtypes", h5("Gene ID type:"),
                     #                                               choices = list("Entrez ID" = "EntrezID", "Gene Symbol" = "GeneSymbol"), selected = "EntrezID")),
                     # #div(style="display: inline-block;vertical-align:top; width: 150px;",selectInput("ddllgra", "Function:",c('mean','median','sd','count','min','max'), selected='mean')),
                     # #div(style="display: inline-block;vertical-align:top; width: 150px;",textInput(inputId="xlimitsmax", label="x-max", value = 0.5))),
            ),
            tabPanel(title = "Enriched Pathways", value = "enrichedPathways",
                     htmlOutput("spacer3"),
                     dataTableOutput("enrichedPathways"),
                     htmlOutput("notes")
            ),
            tabPanel(title = "Gene Hits", value = "signalHits",
                     htmlOutput("spacer4"),
                     tabsetPanel(id = 'signalHits',
                                 #SIGNAL Hits table
                                 tabPanel(title = "SIGNAL Gene Hits", value = "signalHits",
                                          dataTableOutput("signalHits")),
                                 tabPanel(title = "Gene Hits By Iteration", value="geneList",
                                          dataTableOutput("geneList")),
                                 tabPanel(title="Graph: Gene Hits By Iteration", value="geneHitsByIteration",
                                          dataTableOutput("geneHitsTableByIteration"),
                                          plotOutput("geneHitsByIteration")),
                                 tabPanel(title="High Confidence Hits Not Selected by SIGNAL", value="nonSIGNALhits",
                                          dataTableOutput("nonSIGNALhitsTable")),
                                 tabPanel(title = "Pathway Enrichments", value = "pathwayEnrich.cond",
                                          dataTableOutput("pathwayEnrich.cond"))
                     )
            ),
            tabPanel(title = "Network", value = "myNetworkGraph",
                     h4('"Select up to 3 pathways for network graph analysis'), hr(),
                     div(style="display:inline-block",textInput(inputId="mySelection", label="Your selected pathway IDs", value = 0.0)),
                     div(style="display:inline-block",uiOutput("submitGraph")),
                     div(style="display:inline-block",uiOutput("link2Graph")),
                     dataTableOutput("myNetworkGraph")
            ),
            tabPanel(title = "Network Graph", value = "graphViews",
                htmlOutput("spacer5"),
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
                      tabPanel(title = "Network Table", value = "PathNetTable",
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
                     htmlOutput("spacer7"),
                     htmlOutput("downloadFiles"),
                     downloadButton('downloadButton', 'Download all files')
            ),
            tabPanel(title = "User Guide", value = "userGuide",
                     htmlOutput("spacer8"),
                     navlistPanel(id = 'userGuideTab', widths = c(3, 6),
                                  tabPanel(title = "1. Sample dataset", value = "sampleDataset",
                                           uiOutput("sampledataset_page")),
                                  tabPanel(title = "2. Preparing your data", value = "dataPrep",
                                           uiOutput("dataPrep_page")),
                                  tabPanel(title = "3. Running an analysis", value = "runAnalysis",
                                           uiOutput("runAnalysis_page")),
                                  tabPanel(title = "4. Reading SIGNAL results", value = "readResults",
                                           uiOutput("readResults_page")),
                                  tabPanel(title = "5. Saving and securing your analysis", value = "saveAnalysis",
                                           uiOutput("saveAnalysis_page")),
                                  tabPanel(title = "Appendix A: How to segment data into confidence tiers", value = "segmentData",
                                           uiOutput("segmentData_page")),
                                  tabPanel(title = "Appendix B: Data security", value = "security",
                                           uiOutput("dataSecurity_page")),
                                  tabPanel(title = "Appendix C: Running SIGNAL with alternative or bespoke databases", value = "SIGNALcode",
                                           uiOutput("standaloneR"))
                                  )
                     ),
            
            tabPanel(title = "Help", value = "helpUs",
              htmlOutput("spacer9"),
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
      output$spacer8 <- renderUI({
        HTML("<BR>")
      }) 
      output$spacer9 <- renderUI({
        HTML("<BR>")
      }) 
      output$notes <- renderUI({
        HTML("<BR>* <font color=blue>Input High Confidence Hits</font> <BR>* <font color=red>Input Medium Confidence Hits</font><BR>")
      }) 
      # Global environmental variables
      envs <- Sys.getenv()
      env_names <- names(envs)
      
      ## Set up dataDir
      if('SHINY_SERVER_VERSION' %in% env_names){
        dataDir <- '/srv/shiny-server/data/'
      }else{
        dataDir <<- "~/SIGNAL/app/data/"
      }
      
      
      ## Set up downloadDir
      if('SHINY_SERVER_VERSION' %in% env_names){
        outputDir <- '/srv/shiny-server/inputOutputs/SIGNALoutputFiles/'
      }else{
        outputDir <- "~/SIGNAL/app/inputOutputs/SIGNALoutputFiles/"
      }
      userDir <- format(Sys.time(),"%Y%m%d%H%M%S%ms")
      outDir <<- paste0(outputDir, userDir)
      dir.create(outDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
      print(paste0(outputDir, " created"))
      downloadDir <<- paste0(outDir, "/", "SIGNALfilesToDownload")
      currentDir <- getwd()
      dir.create(downloadDir)
      print(paste0(downloadDir, " created"))
      setwd(downloadDir) 
      
      # Read in the input fie
      output$contents <- renderDataTable({
        inFile <- input$file1

        if (is.null(inFile))
          return(NULL)

        # Allow CSV, TXT, EXCEL files ONLY
        if(file_ext(inFile$name) == "csv"){
          data <- read.csv(inFile$datapath, stringsAsFactors = FALSE, header=TRUE)
        }
        else if(file_ext(inFile$name) == "txt"){
          data <- read.table(inFile$datapath, sep = "\t", header=TRUE, fill = TRUE)
        }
        else if(file_ext(inFile$name) == "xlsx"){
          data <- read_excel(inFile$datapath, col_names = TRUE)
        }
        else{
          showModal(modalDialog(title="Input File Type Errors:", HTML("<h3><font color=red>Input file type not supported!<br>Only .CSV, .TXT, .XLSX file allowed</font><h3>")))
          return(NULL)
        }

        # Check for duplicated GeneSymbols
        # if(anyDuplicated(data$GeneSymbol)){
        #   showModal(modalDialog(title="User Input Errors", HTML("<h4><font color=red>Duplicated GeneSymbols were found! <br><br>Please remove the duplicates in the GeneSymobl column and then reload your input file. <br><br>Alternatively, remove the GeneSymbol column from your input file and SIGNAL will map the EntrezID to unique GeneSymbols.</font><h4>", data[duplicated(data$GeneSymbol),]$GeneSymbol)))
        # }
#print(anyDuplicated(data$EntrezID))
        # # Check for duplicated EntrezID
        # if(anyDuplicated(data$EntrezID)){
        #   showModal(modalDialog(title="User Input Errors", HTML("<h4><font color=red>Duplicated EntrezIDs were found! <br><br>Please remove the duplicates and reload your input file.</font><h4>")))
        # }
#print(data$GeneSymbol)       
        ## Complete list of all protein-encoding genes in human genome
        ## ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_types/gene_with_protein_product.txt
        humanGenes <- read.table(file=paste0(dataDir, "HGNC_genes_with_protein_product_EntrezID_geneSymbole_lookup.txt"), sep="\t", header=TRUE)
        
        ## Complete list of all genes in mouse genome
        ## http://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt
        mouseGenes <- read.table(file=paste0(dataDir, "Mouse_genes_with_protein_product_EntrezID_geneSymbol_lookup.txt"), sep="\t", header=TRUE)
        
        ## The two files used above can be updated/regenerated by runing the Rscript 
        ## 'Generate_human_and_mouse_protein_codig_genes_for_SIGNAL.R' in the Rscripts folder
        
        # Check to see if eitehr a 'EntrezID' or a 'GeneSymbol' column is in the input file
        if(!("EntrezID" %in% colnames(data)) && !("GeneSymbol" %in% colnames(data))){
          # Input data do not have EntrezID AND GeneSymbol columns
          showModal(modalDialog(
            title=HTML("<h3><font color=#ff0000>Input file format error!</font></h3>"),
            HTML("Your input file does not contain a required column named 'EntrezID' or 'GeneSymbol'. <br>Please fix your input file and try again!"),
            easyClose = FALSE
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
  
          ############################################################
          # GeneSymbols not recognized (based on GeneBank Entrez gene)
          Unmapped_genes <<- setdiff(data$GeneSymbol, mapped_genes)
        
        
          
          # Save the map and unmapped genes to files into the download folder

          
          data_unmapped <<- data[which(data$GeneSymbol %in% Unmapped_genes), c("GeneSymbol", setdiff(colnames(data), c("GeneSymbol")))]
          setwd(downloadDir)
          fwrite(data_unmapped, file = "unmapped_rows.csv")
          
          # # Save the mapped and unmaped genes for usr to download
          # fwrite(list(mapped_genes), file = "Mapped_genes.csv")
          # fwrite(list(Unmapped_genes), file = "Unmapped_genes.csv")
          setwd(currentDir)
          
          # If no overlapping genes found, catch and handle the error
          if(length(overlappingGenes) == 0){
            showModal(modalDialog(
              title=HTML("<h3><font color=#ff0000>Organism - Gene Set MISMATCH!</font></h3>"),
              #HTML("The organism you selected and the organism from which the input data generated do not match. Please select the correct organism or a different input file, then try again"),
              easyClose = FALSE
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
          
          
          ###Create a modal function for showing a download button:
          myModal <- function() {
            div(id = "test",
                modalDialog(title="Warning:", HTML("<h3><font color=red>Only"), numGeneWithEntrezID,HTML("/"),numGeneInInput, HTML("GeneSymbols have mapped EntrezIDs and will be used in this analysis!</font><h3><br>"),
                            #br(),
                            HTML("<h3><font color=black>You can continue the analysis without the unampped IDs or update your GeneSymbols in the upload file and restart the analysis.</font><h3><br>"),
                            #br(),
                            HTML("<h3><font color=blue>A table of the GeneSymbols that have failed to map to an EntrezID can be viewed in the 'Unmapped Genes' tab or downloaded by clicking on the button below. </font><h3><br>"),
                            #br(),
                            #br(),
                            downloadButton("unmapped_download","Download Unmapped Gene IDs"))
            )
          }
          
          
          

          # Display a warning if one or more input genes have no matching EntrezID due obsolete GeneSymbol
          if ((numGeneWithEntrezID / numGeneInInput) > 0.5 & (numGeneWithEntrezID / numGeneInInput) != 1){
            showModal(myModal())  # Calls the myModal function for a warning message with a download button
            # showModal(modalDialog(title="Warning:", HTML("<h3><font color=red>Only"), numGeneWithEntrezID,HTML("/"),numGeneInInput, HTML("GeneSymbols have mapped EntrezIDs and will be used in this analysis!</font><h3><br>"),
            #                       HTML("Please update your GeneSymbols to match the official symbols in <a href='https://www.ncbi.nlm.nih.gov/gene' target=_blank>Entrez Gene</a> if you want to include ALL in this analysis. Otherwise the following genes will not be used! <br>")))
          }
          else if ((numGeneInInput - numGeneWithEntrezID) > 1){
            showModal(modalDialog(title="Warning:", HTML("<h3><font color=red>Only"), numGeneWithEntrezID,HTML("/"),numGeneInInput, HTML("GeneSymbols have mapped EntrezIDs and will be used in this analysis!</font><h3><br>"),
                                  HTML("Please check if you have the right organism selected for your analysis.")))
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
              easyClose = FALSE
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
        
        updateTabsetPanel(session, "inTabset", selected = "contents")

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
        input.file.columns <<- colnames(siRNA.Score) # Names of input file before SIGNAL analysis

        data <- data %>%
          dplyr::select(GeneSymbol, everything())

        # display the input file dimension
        datatable(data, rownames = FALSE, options = list(paging=TRUE))
      })
      
      output$contents2 <- renderDataTable({
        
        inFile <- input$file1
        
        if (is.null(inFile))
          return(NULL)
        
        # Allow CSV, TXT, EXCEL files ONLY
        if(file_ext(inFile$name) == "csv"){
          data <- read.csv(inFile$datapath, stringsAsFactors = FALSE, header=TRUE)
        }
        else if(file_ext(inFile$name) == "txt"){
          data <- read.table(inFile$datapath, sep = "\t", header=TRUE, fill = TRUE)
        }
        else if(file_ext(inFile$name) == "xlsx"){
          data <- read_excel(inFile$datapath, col_names = TRUE)
        }
        else{
          return(NULL)
        }
        
        
        
        
        if(("EntrezID" %in% colnames(data)) && ("GeneSymbol" %in% colnames(data))) {
          datatable(data.table("The upload file includes both an 'EntrezID' and 'GeneSymbol' column. The analysis will be using the EntrezID ~ GeneSymbol mapping provided by the user."), rownames = FALSE, colnames = NULL, options = list(paging = FALSE, searching = FALSE))
        } else  if(!("EntrezID" %in% colnames(data)) & ("GeneSymbol" %in% colnames(data))){
          if (length(Unmapped_genes) == 0) {
            datatable(data.table("All IDs in your dataset have been succesfully mapped."), rownames = FALSE, colnames = NULL, options = list(paging = FALSE, searching = FALSE))
            } else {
              datatable(data_unmapped, extensions = c('Scroller'), options = list(
                deferRender = TRUE,
                scrollY = 580,
                scroller = TRUE
                )
                )
              
              
                        
              } }
          else if (!("GeneSymbol" %in% colnames(data)) & ("EntrezID" %in% colnames(data))) {
          datatable(data.table("No unmapped genes: EntrezIDs in the upload file with no mapped GeneSymbol were assigned their EntrezID as GeneSymbol. All EntrezIDs will be included in the analysis."), rownames = FALSE, colnames = NULL, options = list(paging = FALSE, searching = FALSE))
          }
        }
        )
      

      
      output$cutoffTypes <- renderUI({
        inFile2 <- input$file1
        if (is.null(inFile2))
          return(NULL)

        # Check validity of the input file
        mtry <- try(read.csv(inFile2$datapath, stringsAsFactors = FALSE, header = TRUE), 
                    silent = TRUE)
        
        if (class(mtry) != "try-error") {
            if(file_ext(inFile2$name) == "csv"){
              data2 <- read.csv(inFile2$datapath, stringsAsFactors = FALSE, header=TRUE)
            }
            else if(file_ext(inFile2$name) == "txt"){
              data2 <- read.table(inFile2$datapath, sep = "\t", header=TRUE, fill = TRUE)
            }
            else if(file_ext(inFile2$name) == "xlsx"){
              data2 <- read_excel(inFile2$datapath, col_names = TRUE)
            }
          } else {
          #showModal(modalDialog(title="User Input Errors:", HTML("<h3><font color=red>Invalid input! Please check your input file.</font><h3>")))
          showModal(modalDialog(title="Input File Type Errors:", HTML("<h3><font color=red>Input file type not supported!<br>Only .CSV, .TXT, .XLSX file format accepted!</font><h3>")))
          
          return(NULL)
        }        
        
        pulldown_types <- c("", colnames(data2))

        selectInput("cutoff_type", "Cutoff Type", pulldown_types)
      })
      
      # Switch to input tab when file is uploaded
      observeEvent(input$file1, {
        inFile2 <- input$file1
        if(is.null(inFile2)) return (NULL)
        updateTabsetPanel(session, "inTabset", selected = "contents")
      })

      # reloads the app
      observeEvent(input$refresh, { 
        session$reload()
      })
      #Secondary cutoff column
      output$secondaryColumn_choice <- renderUI({
        inFile2 <- input$file1
        if (is.null(inFile2) | is.null(input$cutoff_type) |   input$secondaryCutoff == FALSE)
          return(NULL)
        
        # if(input$cutoff_type == "")
        #   return(NULL)
        
        
        
        data2 <- read.csv(inFile2$datapath, stringsAsFactors = FALSE,header = TRUE)
        

        pulldown_types <- c("", colnames(data2))
        
        pulldown_types_2 <- setdiff(pulldown_types, input$cutoff_type)
        
        selectInput("secondaryColumn_menu", "Column to Use for Secondary Criteria", pulldown_types_2)
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
          scriptDir <- "~/SIGNAL/app/Rscripts/"
        }

        if('SHINY_SERVER_VERSION' %in% env_names){
          inputDir <- '/srv/shiny-server/inputOutputs/SIGNALinputFiles/'
        }else{
          inputDir <- "~/SIGNAL/app/inputOutputs/SIGNALinputFiles/"
        }

        if('SHINY_SERVER_VERSION' %in% env_names){
          outputDir <- '/srv/shiny-server/inputOutputs/SIGNALoutputFiles/'
        }else{
          outputDir <- "~/SIGNAL/app/inputOutputs/SIGNALoutputFiles/"
        }

        # # To keep a copy of html files for iframe to access
        if('SHINY_SERVER_VERSION' %in% env_names){
          wwwDir <<- '/srv/shiny-server/www/'
        }else{
          wwwDir <<- "~/SIGNAL/app/www/"
        }
        
        # Get organism name from user input
        organism <<- input$organism
        organismAbbr <- ifelse(grepl("human", tolower(organism)), 'hsa', 'mmu')

        ## Source other codes depending on whether this is used
        # as a standalone tool or on AWS webservice
        if('SHINY_SERVER_VERSION' %in% env_names){
          source("/srv/shiny-server/Rscripts/pathway_iteration.R", local = TRUE)
        }else{
          source("~/SIGNAL/app/Rscripts/pathway_iteration.R", local = TRUE)
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
        
        dir_begin <<- ifelse('SHINY_SERVER_VERSION' %in% env_names, paste0("/srv/shiny-server/data/Networks/String."), paste0("~/SIGNAL/app/data/Networks/String."))
        
        if(input_netowrk == "STRING: Experimental & Database"){
          
          # loads igraph for experimental and database source selections
          load(paste0(dir_begin, tolower(organism), ".exp_and_data", conf.level, ".Rdata"))
          G <- get(paste0("String.", tolower(organism), ".exp_and_data", conf.level))
          
        }
        else if(input_netowrk == "Advanced Options"){
          files2load <<- paste0(dir_begin, tolower(organism), ".", network_InteractionSources, conf.level, ".Rdata")
          vars2load <<- paste0("String.", tolower(organism), ".", network_InteractionSources, conf.level)
          for (i in 1:length(files2load)) assign(gsub(".*/","",files2load[i]), load(files2load[i]))
          all.G = list()
          for(i in 1:length(vars2load)) all.G[[i]] = get(vars2load[i])
          empty_func <- function(x) length(unlist(x)) == 0L
          empty.G.check = sapply(all.G, edge.attributes) %>% lapply(`[`, 'datasource') %>%
            lapply(empty_func) %>% unlist() %>% which()
          count.G.check = which(sapply(all.G, vcount) <= 0)
          all.checks = union(count.G.check, empty.G.check)
          
          datasources = sapply(strsplit(vars2load, '[,.]'), '[[', 3)
          
          # check if the selected was experimental and database or only one of them
          if(length(match(datasources, c("experimental", "database"))) == 2){
            load(paste0(dir_begin, tolower(organism), ".exp_and_data", conf.level, ".Rdata"))
            G <- get(paste0("String.", tolower(organism), ".exp_and_data", conf.level))
          }
          else{
            if(length(all.checks)>0){
              all.G = all.G[-all.checks]
              datasources = datasources[-all.checks]
            }
            if(length(all.G)==0){
              showModal(modalDialog(title="Warning:", HTML("<h3><font color=red>Criteria produced empty network. Session will restart.</font><h3>"),
                                    easyClose = FALSE))
              Sys.sleep(3)
              session$reload()
            } else if(length(all.G)==1){
              G <- all.G[[1]]
            } else{
              
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
        }
        #Selected_STRINGnetwork.igraph <- G
        message("Networks Loaded")

        #message(networkType)
        use.only.commnected.components <- c('Yes')

        # Remove all files in the outputDir
        unlink(outputDir, recursive = FALSE)
        
        # Create user-specific directory using system time
        # userDir <- format(Sys.time(),"%Y%m%d%H%M%S%ms")
        # outDir <<- paste0(outputDir, userDir)
        # dir.create(outDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
        setwd(outDir)
        
        # Set the output file name
        inputFile <- input$file1
        inputFileName <<- inputFile$name
        inputFilePrefix <<- tools::file_path_sans_ext(inputFileName)

        outputFileName <<- paste0(inputFilePrefix, "_", "SIGNALoutput_ALL.csv")

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
        
        pathwayData <<- read.csv(file = paste0(dataDir, "Pathways/KEGG", KEGG.version, "_", organism, "_", pathway.type, ".csv"))

        proxyScore <- input$cutoff_type
        
        
        #########################################
        ########## Create data tier of gene lists
        #########################################
        cutoffHigh <<- as.numeric(input$cutoff_valueH)
        cutoffMed <<- as.numeric(input$cutoff_valueM)
        cutoffType <<- input$cutoff_type
        
        secondary_column <<- input$secondaryColumn_menu
        secondary_direction <<- input$secondaryDirection
        secondary_value <<- as.numeric(input$secondaryValue)
        
        

        
        
        #Ensure cutoff column is as numeric
        siRNA.Score <- 
        
        #Depending on the differnce between the high conf cutoff and the mid conf cutoff assign criteria based on "greater than" or "less than"
        if((cutoffHigh - cutoffMed) > 0){
          siRNA.Score <- siRNA.Score %>%
            mutate(InputCategory = ifelse(get(cutoffType, as.environment(siRNA.Score)) >= cutoffHigh, "HighConf",
                                               ifelse(get(cutoffType, as.environment(siRNA.Score)) >= cutoffMed & get(cutoffType, as.environment(siRNA.Score)) < cutoffHigh,"MedConf",
                                                      "LowConf")))
        }else{
          siRNA.Score <- siRNA.Score %>%
            mutate(InputCategory = ifelse(get(cutoffType, as.environment(siRNA.Score)) <= cutoffHigh, "HighConf",
                                               ifelse(get(cutoffType, as.environment(siRNA.Score)) <= cutoffMed & get(cutoffType, as.environment(siRNA.Score)) > cutoffHigh,"MedConf",
                                                      "LowConf")))
        }
        
        #### Depnding on if there's another cirteria check that all hits meet criteria
        additonalCriteria <- input$secondaryCutoff
        
        
        if(additonalCriteria){
          if (secondary_direction == "\u2265"){
            siRNA.Score <- siRNA.Score %>%    #greater than
              mutate(InputCategory = ifelse(get(secondary_column, as.environment(siRNA.Score)) >= secondary_value, InputCategory, "LowConf"))
          } else {
            siRNA.Score <- siRNA.Score %>%   #less than
              mutate(InputCategory = ifelse(get(secondary_column, as.environment(siRNA.Score)) <= secondary_value, InputCategory, "LowConf"))
          }
        }
        
        ## To include background genes 
        includeBackground <- input$includeBackground
        
        if(includeBackground){
          df_background <- data.frame(EntrezID = df_backgroundGenes$EntrezID, GeneSymbol=df_backgroundGenes$GeneSymbol, InputCategory = rep("Background", nrow(df_backgroundGenes)))
          # combined input data and background data
          myList <- list(siRNA.Score, df_background)
          siRNA.Score <- rbindlist(myList, fill = TRUE)
        }
        


        
        # To assign genes that never get filtered out
        if(input$manualHits != "") {
          manual.hits <- unlist(str_split(str_replace_all(str_replace_all(input$manualHits, "[\r\n]" , ","), "[\\s]" , ""), ","))
          

          
          if(input$manualIDtypes == "EntrezID"){
            siRNA.Score <- siRNA.Score %>%   
              mutate(InputCategory = ifelse(EntrezID %in% manual.hits, "AssignedHit", InputCategory))
          } else {
            siRNA.Score <- siRNA.Score %>%  
              mutate(InputCategory = ifelse(tolower(GeneSymbol) %in% tolower(manual.hits), "AssignedHit", InputCategory))
          }
        }
        

        
        ###### For analysis

        cutoffType <<- input$cutoff_type
        proxyScore <- "InputCategory"
        iteration <- 1
        counter <- TRUE

        # Get a copy of the original list of high-confidence genes
        originalHits <- siRNA.Score$GeneSymbol[siRNA.Score$InputCategory %in% c("HighConf", "AssignedHit")]

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

        # Perform iterative SIGNAL analysis
        
        while (counter) {
          
          ## Show progress
          # Increment the progress bar, and update the detail text.
          progress$inc(1/(iteration*3))
          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(0.1)
          
          Hits <- siRNA.Score$EntrezID[siRNA.Score[[proxyScore]] == "HighConf" | siRNA.Score[[proxyScore]] == "AssignedHit"]

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
          siRNA.Score[[kName1]][siRNA.Score$InputCategory == "AssignedHit"] <- "HighConf"  #Mnaually assigned hits
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
          siRNA.Score[[nName2]][siRNA.Score$EntrezID %in% gNames2 & siRNA.Score[[kName1]] %in% c("MedConf", "HighConf", "AssignedHit")] <- "HighConf"
          
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

        ## Save the results into output files in the SIGNALoutputFiles folder
        signal.Out <<- merge(siRNA.Score, pathDF, by = "GeneSymbol", all = T)

        message(getwd(), "#####")

        write.csv(signal.Out, file = outputFileName, row.names = F)
        final_enriched_pathway_file <- paste0(pathway.type, "_SIGNAL_enrichment_final", ".csv")
        write.csv(pathEnrich, file = final_enriched_pathway_file, row.names = F)
        
        # Pass the pathway list to build networkGraph
        sigPathways <<- pathEnrich[,c("Pathway", "Genes", "HitGenes")]

        completed <- TRUE
     
      ###################################
      ### Generate Condensed Output Files
      ###################################
      
      SIGNALoutput <- read.csv(outputFileName, stringsAsFactors = F)
        cutoffType <<- input$cutoff_type
        cutoffHigh <<- as.numeric(input$cutoff_valueH)
        cutoffMed <<- as.numeric(input$cutoff_valueM)
        
        #################
        # Set High and Low confidence Label and SIGNAL Label
        ##################
        
        FinalIterationNetworkColumn <- paste0("Network.class.iteration", iterationNum)
        
        SIGNALoutput <- SIGNALoutput %>%
          mutate(SIGNALhit = ifelse(get(FinalIterationNetworkColumn, envir = as.environment(SIGNALoutput)) == "HighConf", 
                                    "Yes",
                                    ""))
        

        ################
        # Generate Matrices for SIGNAL Hits, high conf hits, med conf hits
        #Get filtered SIGNALhits                                           # This is where I put together what is considered a "hit" by SIGNAL (IAM). Any gene that had a score of 1 in the last network analysis step
        
        SIGNALhits <<- filter(SIGNALoutput, SIGNALhit == "Yes")
        SIGNALhits.matrix <- matrix(SIGNALhits$EntrezID)                      # Created a matrix of all the genes that are "hits" 
        
        # Get filtered SIGNAL High/Med Conf
        SIGNALhits.highConf <- filter(SIGNALoutput, SIGNALhit == "Yes" 
                                      & InputCategory %in% c("HighConf", "AssignedHit"))
        SIGNALhits.highConf.matrix <- matrix(SIGNALhits.highConf$EntrezID)
        SIGNALhits.highConf.matrix.GS <- matrix(SIGNALhits.highConf$GeneSymbol)
        
        SIGNALhits.medConf <- filter(SIGNALoutput, SIGNALhit == "Yes" 
                                     & InputCategory == "MedConf")
        SIGNALhits.medConf.matrix <- matrix(SIGNALhits.medConf$EntrezID)
        SIGNALhits.medConf.matrix.GS <- matrix(SIGNALhits.medConf$GeneSymbol)
        
        #numTotal variabale to be used in GeneList tab to generate graph
        
        numTotal  <<- length(SIGNALhits.highConf.matrix) + length(SIGNALhits.medConf.matrix)
        
        #############################################################     # Using the methods from CARD (only swithced to EntrezID over GeneSymbol)
        #            PageRank Algorithm to All Genes
        #############################################################               
        # AVERAGE DUPLICATED ROWS
        siRNA.Score <- SIGNALhits[which(!is.na(SIGNALhits$EntrezID)),]    
        OverallDegree <- degree(G)
        
        Screen_Genes.for.network.analysis <- intersect(siRNA.Score$EntrezID, V(G)$name)
        
        Graph <- induced.subgraph(G, Screen_Genes.for.network.analysis)
        
        #Screen_Genes.for.network.analysis <- intersect(siRNA.Score$EntrezID,V(G)$name)
        
        #############################################################
        #            Select sub-graphs from hit genes
        #############################################################
        SubGraph <<- Graph
        
        if(length(E(SubGraph))==0 | length(V(SubGraph))==0){
          showModal(modalDialog(title="Warning:", HTML("<h3><font color=red>Subgraph produced empty network due to low number of hit genes. Session will restart.</font><h3>"),
                                easyClose = FALSE))
          Sys.sleep(3)
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
        SIGNALhits.merge = SIGNALhits[, c("EntrezID", "GeneSymbol", "InputCategory", "Pathway", "SIGNALhit")]
        GraphNodesHit <-  merge(GraphNodesHit, SIGNALhits.merge,                                                  
                               by.x = "EntrezID", by.y = "EntrezID", all.x = T)
        
        ################    
        # Now getting a data frame for the edges and a data frame for the nodes
        ################ EDGE INFO and NODE INFO are sent to GLOBAL ENVIRONMENT to be used by RANKING
        
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
        #Merge SIGNALhits to NodeInfo
        Scores_and_nodes <- merge(NodeInfo[, c("GeneMappingID", "GeneSymbol")],                      #Pairing up the "Node Info" with the gene info (such as gene symbol and groupings)
                                  SIGNALhits, 
                                  by.x = "GeneSymbol", by.y = "GeneSymbol", all.y = T)
        ###################
        ######Create Data frame with Pathway interactome
        #Convert Target values to pathways

        EdgeInfo.TargetPathways <- merge(EdgeInfo, Scores_and_nodes[, c("GeneMappingID", "Pathway")],
                                         by.x = "target", by.y = "GeneMappingID", all.x = T)
      
        #EdgeInfo.TargetPathways <- merge(EdgeInfo, Scores_and_nodes[, c("GeneMappingID", "Pathway")],
        #                                 # by.x = "target", by.y = "GeneMappingID",  allow.cartesian=TRUE)
        
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
        #Align column names and stack data frames                                                         
        #The analysis assumed directionality of interactions, but we're ignoring it here, so combining the "target" and Source" to one dataframe.
        colnames(EdgeInfo.SourcePathways_Sum) <- c("source", "Pathway")
        EdgePathways_stacked <- rbind(EdgeInfo.TargetPathways_Sum, EdgeInfo.SourcePathways_Sum)
        
        
        #Aggregate
        EdgePathways_stacked_Sum <<- EdgePathways_stacked %>% 
          group_by(source) %>% summarise(Pathway = toString(Pathway))
        
        ##Add counts to pathways (number of genes in list that are part of each pathway)
        # EdgePathways_stacked_Sum$Pathway.counts <- NA
        # 
        # for (i in 1:length(EdgePathways_stacked_Sum$Pathway)) {
        #   temp.string <<- unlist(strsplit(EdgePathways_stacked_Sum$Pathway[i], ", "))
        #   for (j in 1:length(temp.string)) {
        #     name.j <<- temp.string[j]
        #     count.j <<- length(grep(name.j, temp.string, fixed = T))
        #     out <- paste0(name.j, " (", count.j, ")")
        #     EdgePathways_stacked_Sum$Pathway.counts[i] <-  ifelse(j == 1, out, paste(out, EdgePathways_stacked_Sum$Pathway.counts[i], sep = ", ")) 
        #   }
        # }
        # 
        # 
        # 
        # 
        # ########Remove Duplicates of Pathway names
        # 
        # EdgePathways_stacked_Sum$NetworkGenePathways <- sapply(strsplit(EdgePathways_stacked_Sum$Pathway.counts, ", ", fixed = TRUE), function(x) 
        #   paste(unique(x), collapse = ", "))
        # 
        # ###########Order from pathway with highest number of genes to lowest
        # for (k in 1:length(EdgePathways_stacked_Sum$Pathway)) {
        #   pathway.string <- unlist(strsplit(EdgePathways_stacked_Sum$NetworkGenePathways[k], ", "))
        #   pathway.values <- as.numeric(gsub("[\\(\\)]", "", regmatches(pathway.string, gregexpr("\\(.*?\\)", pathway.string))))
        #   names(pathway.values) <- order(pathway.values, decreasing = T)
        #   EdgePathways_stacked_Sum$NetworkGenePathways[k] <- paste(pathway.string[as.numeric(names(pathway.values))], collapse = ", ")
        # }
        
        writing.pathways <- function(edge.df){
          edge.df$Pathway.counts <- NA
          L = sapply(edge.df$Pathway, strsplit, ", ")
          names(L) = edge.df$source
          tab = lapply(L, table)
          tab = lapply(tab, sort, decreasing=T)
          
          l.pathway.counts = lapply(tab, function(x){
            pathway.counts = c()
            for(i in 1:length(x)){
              x.name = names(x)[i]
              pathway.counts = append(pathway.counts, paste0(x.name, " (", as.numeric(x[i]), ")"))
            }
            return(pathway.counts)
          })
          
          l.pathway.cond = lapply(l.pathway.counts, paste0, collapse=', ')
          source.name = names(l.pathway.cond)

          network.paths = as.vector(unlist(l.pathway.cond))
          
          edge.df = data.frame("source" = source.name, "NetworkGenePathways" = network.paths)
          return(edge.df)
        }
        
        EdgePathways_stacked_SumandValue = writing.pathways(EdgePathways_stacked_Sum)
        
        
        #######Add in entrezID for source genemappings
        EdgePathways_stacked_EntrezID <- merge(Scores_and_nodes[, c("GeneMappingID", "EntrezID")], EdgePathways_stacked_SumandValue,
                                               by.x = "GeneMappingID", by.y = "source",
                                               all.y = T)
        
        #######Add to SIGNAL output file
        SIGNALoutput <- merge(SIGNALoutput, EdgePathways_stacked_EntrezID[, c("EntrezID", "NetworkGenePathways")],
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
        
        #Update Names                                                                                     #Now a dataframe is being created where each gene selected by SIGNAL (or part of the highlighted groups) has a list "Ntwrk.all" that lists all other genes from TRAIGE that it is predicted to interact with.
        colnames(Edge_summary) <- c("GeneSymbol", "InteractingGenes")
        
        #Merge with scores 
        SIGNALoutput <- merge(SIGNALoutput, Edge_summary, by.x = "GeneSymbol", by.y = "GeneSymbol", all.x = T)
        
        # table where high confidence genes are not in the final iteration of the SIGNAL analysis
        non.signal.dat <- filter(SIGNALoutput, InputCategory=="HighConf" & SIGNALhit=='')
     
        if(dim(non.signal.dat)[1]){
          non.signal.dat$SIGNALhit = "No"
          # non.signal.dat <- non.signal.dat[,c("EntrezID", "GeneSymbol", "InputCategory", "SIGNALhit", "Pathway", "InteractingGenes", "NetworkGenePathways")]
          #non.signal.dat <- non.signal.dat[,c("EntrezID", "GeneSymbol", "InputCategory", "SIGNALhit")]
          non.signal.dat <- non.signal.dat[ ,c(unique(c("EntrezID", "GeneSymbol", input.file.columns, "InputCategory",  "SIGNALhit")))]  ## remove added on columns
          
          # output non.signal.dat table
          output$nonSIGNALhitsTable <- renderDataTable({
            return(non.signal.dat)
          })
        }
        
        ####################
        ### Create Condensed Output File
        SIGNALoutput.condensed <<- SIGNALoutput[SIGNALoutput$SIGNALhit == "Yes", c("EntrezID", "GeneSymbol", "InputCategory", "SIGNALhit", "Pathway", "InteractingGenes", "NetworkGenePathways")]

        ## Change the column name from 'HitGenes' to 'SIGNALhits'
        names(SIGNALoutput.condensed)[names(SIGNALoutput.condensed) == "SIGNALhit"] <- "SIGNALhits"
        
        
        ########################
        ######## Pathway Output
        #########################
        FinalEnrichment.df <- pathEnrich
        
        
        #Genrate columns with high confidence and med confidence (based on input) genes of each pathway.
        FinalEnrichment.df$HighScoreHitGenes <- NA
        FinalEnrichment.df$HighScoreHitGenesNames <- NA
        FinalEnrichment.df$MedScoreGenesNames <- NA
        
        for (i in 1:length(FinalEnrichment.df$Genes)) {
          temp.path.string <- unlist(strsplit(FinalEnrichment.df$HitGeneNames[i], ", "))
          out.HC <- paste(intersect(temp.path.string, SIGNALhits.highConf.matrix.GS),collapse = ", ")
          out.HC.count <- length(intersect(temp.path.string, SIGNALhits.highConf.matrix.GS))
          out.MC <- paste(intersect(temp.path.string, SIGNALhits.medConf.matrix.GS),collapse = ", ")
          FinalEnrichment.df$HighScoreHitGenes[i] <- out.HC.count
          FinalEnrichment.df$HighScoreHitGenesNames[i] <- out.HC
          FinalEnrichment.df$MedScoreGenesNames[i] <- out.MC
        }
        
        
        ########### Generate Enrichment Score
        FinalEnrichment.df$EnrichScore <- NA
        
        for (i in 1:length(FinalEnrichment.df$Pathway)) {
          GeneHitGeneRatio <- FinalEnrichment.df$HitGenes[i] / FinalEnrichment.df$Genes[i]
          HighConfHitGeneRation <-  FinalEnrichment.df$HighScoreHitGenes[i] / FinalEnrichment.df$HitGenes[i]
          FinalEnrichment.df$EnrichScore[i] <- round(((GeneHitGeneRatio + HighConfHitGeneRation) / 2), 3)
        }
        
        FinalEnrichment.condensed <- FinalEnrichment.df[, c("Pathway", "pVal", "pValFDR", "pValBonferroni", "Genes", "HitGenes", "HighScoreHitGenes", "HighScoreHitGenesNames", "MedScoreGenesNames", "EnrichScore")]
        
        ## Change the 'HitGenes' to 'SIGNALhits'
        names(FinalEnrichment.condensed)[names(FinalEnrichment.condensed) == "HitGenes"] <- "SIGNALhits"
        
        ############# Write files to new Directory
        #downloadDir <- paste0(outDir, "/", "SIGNALfilesToDownload")
        #dir.create(downloadDir)
        # setwd(downloadDir)
        setwd(downloadDir)
        
        SIGNAL.cond.output.name <<- paste0(inputFilePrefix, "_", "SIGNALhits.csv")
        Enrichment.cond.output.name <- paste0(inputFilePrefix, "_", "SIGNALenrichment.csv")
        
        fwrite(SIGNALoutput.condensed, file = SIGNAL.cond.output.name)
        fwrite(FinalEnrichment.condensed, file = Enrichment.cond.output.name)
        
        if(dim(non.signal.dat)[1]){
          fwrite(non.signal.dat, file = paste0(inputFilePrefix, "_", "nonSIGNALhits.csv"))
          # write.csv(signal.Out, file = outputFileName, row.names = F)
        }

      ######################
      ## Change 'HitGenes' to 'SIGNALhits' in the display table under 'Enriched Pathways' tab
      names(pathEnrich)[names(pathEnrich) == "HitGenes"] <- "SIGNALhits" 
        
        
        
      ## Switch to 'Enriched Pathways' tab and display partial results
      observe({
        if(completed) {
          updateTabsetPanel(session, "inTabset", selected = "enrichedPathways")

          output$enrichedPathways <- renderDataTable({
            options = list(autoWidth = TRUE, scrollX = TRUE,
                           columnDefs = list(width = "100%", targets = "HitGeneNames"))

            # Used to add hyperlink to KEGG pathway
            fontBlue <- function(val) {
              sprintf('<font color=blue>%s</font>',val)
            }

            fontRed <- function(val) {
              sprintf('<font color=red>%s</font>',val)
            }

            # Add a hyperlink to KEGG pathway
            createLink <- function(val1, val2, val3) {
              sprintf('<a href="https://www.kegg.jp/kegg-bin/mcolor_pathway?%s0%s" target="_blank" class="btn btn-primary">%s</a>', val1, val2, val3)
            }

            # Add a hyperlink to KEGG mapper
            link2KEGGmapper <- function(organismAbbr, pathwayID,  myGeneLabels, pathwayName) {
              # Check the data submitted by the form using http post method
              #sprintf('<form target="_blank" enctype="multipart/form-data" method="post" action="http://localhost/cgi-bin/display_form_data.cgi">
              # Create a form for each datatable row
              sprintf('<form onsubmit="return confirm(\'You are going to a non-NIH website! This external link provides additional information that is consistent with the intended purpose of this site. NIH cannot attest to the accuracy of a non-federal site. Linking to a non-federal site does not constitute an endoresment by NIH or any of its employees of the sponsors or the information and products presented on the site. You will be subject to the destination site privacy policy when you follow the link.\');" target="_blank" enctype="multipart/form-data" method="post" action="https://www.kegg.jp/kegg-bin/mcolor_pathway?" target="_blank">
                       <input type="hidden" name="map" value="%s%s">
                       <input type="hidden" name="unclassified" value="%s">
                       <input type="hidden" name="s_sample" value="color">
                       <input type="hidden" name="mode" value="color">
                       <input type="hidden" name="reference" value="white">
                       <input type="submit" style="font-face: \'Comic Sans MS\'; font-size: larger; color: teal; background-color: powderblue; border: 0 none;"value="%s"></form>', organismAbbr, pathwayID, myGeneLabels, pathwayName)
            }
            

            # Used to add text color to GeneSymbols indicate whether they are in the original hit or identified by SIGNAL
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
                  myRedGeneLabels <- paste(myRedGeneLabels, myRedGeneLabel, sep="\n")
                  myRedGeneIDs <- capture.output(cat(myRedGeneIDs, myRedGeneLabel))
                }
              }

              # Add hyperlink to KEGG pathway database
              pathwayName <- pathEnrich[i,][1]
              pathwayID <- pathwayData$PathwayID[match(pathwayName, pathwayData$PathwayName)]
              pathwayID_full <- ifelse(nchar(pathwayID) == 3, paste0("00", pathwayID), paste0("0", pathwayID))  #Since kegg requires the zeros before the number as well add back zeros that that the total number of characters is 5
              mapperHeader <- capture.output(cat("#", organismAbbr,	"CLP/CMP\tBlast_phase\tAll"))
              myGeneLabels <- paste(mapperHeader, stri_replace_all_fixed(myBlueGeneLabels, " ", ""), "\n", stri_replace_all_fixed(myRedGeneLabels, " ", ""), sep = "")
              pathEnrich[i,][1] <- link2KEGGmapper(organismAbbr, pathwayID_full, myGeneLabels, pathwayName)

              # Display the original hits(BLUE) first, followed by the hits picked up by SIGNAL (RED)
              #myGene <- paste(myBlueGene, myRedGene, sep = "")
              myRedGene <- sub('.', '', myRedGene)
  
              # Remove <br> if no blue genes are associated with the pathway
              myGene <- paste(myBlueGene, myRedGene, sep = "<br>")
              myGene <- sub('^<br>', '<', myGene)
  
              myGene <- substring(myGene, 2)
              pathEnrich[i,7] <- myGene
            }
            #View(pathEnrich)
            # Chang column name from 'Genes' to 'PathwayGenes'
            colnames(pathEnrich)[which(names(pathEnrich) == "Genes")] <- "PathwayGenes"
            return(pathEnrich)
            completed2 <- TRUE
          }, escape = FALSE)
        }
      })

      # Create the 'Gene List' tab
      # output$geneHits <- renderUI({
      #   updateTabsetPanel(session, "inTabset", selected = "geneHits")
        
        output$signalHits <- renderDataTable({
          dat <- datatable(SIGNALoutput.condensed, rownames = FALSE, options = list(paging=T, autoWidth = F, scrollX = F
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
        
        geneHitsToPlot <- reactive({
          EnrichColumns.index <- NULL
          SIGNALrenamed <- SIGNALoutput
          
          
          for (i in 1:iterationNum){
            pathColumn <- which(colnames(SIGNALrenamed) == paste0("KEGG.class.iteration", i))
            netColumn <- which(colnames(SIGNALrenamed) == paste0("Network.class.iteration", i))
            colnames(SIGNALrenamed)[pathColumn] <- paste0("PathwayEnrichment", i)
            colnames(SIGNALrenamed)[netColumn] <- paste0("NetworkEnrichment", i)
            EnrichColumns.index <- c(EnrichColumns.index, pathColumn, netColumn)
          }
          
          SIGNALiterations <- SIGNALrenamed[, c(1:(min(EnrichColumns.index)-2), EnrichColumns.index, (max(EnrichColumns.index)+1):length(SIGNALoutput))]
          Iteration.index <- c(which(colnames(SIGNALiterations) == "InputCategory"), grep('NetworkEnrichment', names(SIGNALiterations)))
          
          
          totalRow <- data.frame(matrix(NA,1,length(SIGNALiterations)))
          colnames(totalRow) <- colnames(SIGNALiterations)
          hitsDataFrame <<- data.frame(matrix(0, length(Iteration.index), 4))
          colnames(hitsDataFrame) <- c('Iteration', 'Total', 'High-conf', 'Med-conf')
          
          for (l in 1:length(Iteration.index)){
            totalHits <- length(which(SIGNALiterations[Iteration.index[l]] == "HighConf"))
            totalHighConf <- length(which(SIGNALiterations[Iteration.index[l]] == "HighConf" & SIGNALiterations$InputCategory == "HighConf"))
            totalMedConf <- length(which(SIGNALiterations[Iteration.index[l]] == "HighConf" & SIGNALiterations$InputCategory == "MedConf"))
            totalRow[1, Iteration.index[l]] <- totalHits
            hitsDataFrame[l, ] <- c(l - 1, totalHits, totalHighConf, totalMedConf)
          }
          return(hitsDataFrame)
        })
        
        SIGNALiterations <- reactive({
          EnrichColumns.index <- NULL
          SIGNALrenamed <- SIGNALoutput
          
          
          for (i in 1:iterationNum){
            pathColumn <- which(colnames(SIGNALrenamed) == paste0("KEGG.class.iteration", i))
            netColumn <- which(colnames(SIGNALrenamed) == paste0("Network.class.iteration", i))
            colnames(SIGNALrenamed)[pathColumn] <- paste0("PathwayEnrichment", i)
            colnames(SIGNALrenamed)[netColumn] <- paste0("NetworkEnrichment", i)
            EnrichColumns.index <- c(EnrichColumns.index, pathColumn, netColumn)
          }
          
          SIGNALiterations <- SIGNALrenamed[, c(1:(min(EnrichColumns.index)-2), EnrichColumns.index, (max(EnrichColumns.index)+1):length(SIGNALoutput))]
          Iteration.index <- c(which(colnames(SIGNALiterations) == "InputCategory"), grep('NetworkEnrichment', names(SIGNALiterations)))
          
          
          totalRow <- data.frame(matrix(NA,1,length(SIGNALiterations)))
          colnames(totalRow) <- colnames(SIGNALiterations)
          hitsDataFrame <<- data.frame(matrix(0, length(Iteration.index), 4))
          colnames(hitsDataFrame) <- c('Iteration', 'Total', 'High-conf', 'Med-conf')
          
          for (l in 1:length(Iteration.index)){
            totalHits <- length(which(SIGNALiterations[Iteration.index[l]] == "HighConf"))
            totalHighConf <- length(which(SIGNALiterations[Iteration.index[l]] == "HighConf" & SIGNALiterations$InputCategory == "HighConf"))
            totalMedConf <- length(which(SIGNALiterations[Iteration.index[l]] == "HighConf" & SIGNALiterations$InputCategory == "MedConf"))
            totalRow[1, Iteration.index[l]] <- totalHits
            hitsDataFrame[l, ] <- c(l - 1, totalHits, totalHighConf, totalMedConf)
          }
          totalRow[1,1] <- "Total"
          SIGNALiterations <- rbind(totalRow, SIGNALiterations)
          return(SIGNALiterations)
        })
        
        output$geneList <- renderDataTable({
            
            
            
            # View the dataframe 
            #geneHitsToPlot <<- data.frame(hitsDataFrame)
            
            # Highlight the 'Total' row using formatStyle()
            dat <- datatable(SIGNALiterations(), rownames = FALSE, options = list(paging=T, autoWidth = F, scrollX = F, 
                                                                                columnDefs = list(list(width = '200px', length = '400px',
                                                                                                         targets = c((length(SIGNALiterations())-2), (length(SIGNALiterations())-1)),
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
            dat1 <- datatable(geneHitsToPlot(), options = list(paging = FALSE, searching = FALSE, rownames = FALSE))
            return (dat1)
          })
          
          # Create plots showing the numbers of gene hits by iteration
          output$geneHitsByIteration <- renderPlot({
            geneHitsToPlot.melt <<- melt(geneHitsToPlot(), id.vars = "Iteration")
            ggplot(data = geneHitsToPlot.melt, aes(x = as.integer(Iteration), y = as.numeric(value), group = variable, color = variable)) +
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
        
        updateTabsetPanel(session, "inTabset", selected = "myNetworkGraph")

        colnames(sigPathways) <- c("Pathway", "PathwayGenes", "SIGNALhits")

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
            message("Inside NetworkGraph2")
            
            source(paste0(scriptDir, "config_jsons.R"), local = TRUE)
            #source(paste0(scriptDir, "Ranking_plusComments_v3.R"), local = TRUE)
            source(paste0(scriptDir, "Ranking_source.R"), local = TRUE)
            progress1$inc(1/2)
            
            # Need to catch error to allow reload the app
            
            Generate_NetworkGraph(selectedRows, organism, G)
            
            # Writing fully generated network files for download
            # fwrite(rbindlist(json_2df), file = paste0(inputFilePrefix, "_", "Second_Degree_Network.csv"))
            # fwrite(rbindlist(json_1df), file = paste0(inputFilePrefix, "_", "First_Degree_Network.csv"))
            
            
            print(paste0("Network Generation Time: ", Sys.time() - startB))
            
            # reactive statement for listing output files
            # downloadOut <- reactive({
            #   outputFiles = list.files(path = './')
            #   out <- c("<br><b>All files from your SIGNAL analysis for download:</b><br>")
            # 
            #   if(length(outputFiles) > 0){
            #     for(i in outputFiles){
            #       out <- paste(out, i, sep = "<br>")
            #     }
            #   }
            # 
            #   out <- paste0(out, "<br><br>")
            #   return(out)
            # })
            
            downloadOut <- reactive({
              outputFiles = list.files(path = './')
              # reordering files so that enrichment and SIGNAL hits files come first and second  and as a seperate list on top
              enrich.spot = grep("SIGNALenrichment", outputFiles)
              hits.spot = grep("_SIGNALhits", outputFiles)
              main_outputFiles = c(outputFiles[c(hits.spot, enrich.spot)])
              
              
              #Select files for non hits
              nonhits.spot = grep("_nonSIGNALhits", outputFiles)
              nonhits_outputFiles = c(outputFiles[c(nonhits.spot)])
              
              #Select files for networks
              network.spot = grep("SIGNALnetwork", outputFiles)
              clickpath.spot = grep("Clicked_Pathways", outputFiles)
              network_outputFiles = c(outputFiles[c(network.spot, clickpath.spot)])
              
              #select all other files not selected above
              
              other_outputFiles = c(outputFiles[-c(enrich.spot, hits.spot, nonhits.spot, network.spot, clickpath.spot)])
              
              
              
              out <- c("<br><b>All files from your SIGNAL analysis for download:</b><br> <br><b><font color=blue>SIGNAL selected hits and pathway enrichments:</font></b><br>")
              
              
              #Main output files
              if(length(main_outputFiles) > 0){
                for(i in main_outputFiles){
                  out <- paste(out, paste0("<i>", i, "</i>"), sep = "<br>")
                }
              }
              
              #nonhits output files
              if(length(nonhits_outputFiles) > 0){
                out <- paste(out, c("<br><b><font color=blue>High confidence hits that were not selected by SIGNAL:</font></b><br>"), sep = "<br><br>")
                
                for(i in nonhits_outputFiles){
                  out <- paste(out, paste0("<i>", i, "</i>"), sep = "<br>")
                }
              }
              
              #nnetwork output files
              if(length(network_outputFiles) > 0){
                out <- paste(out, c("<br><b><font color=blue>Selected networks and clicked pathway files:</font></b><br>"), sep = "<br><br>")
                
                for(i in network_outputFiles){
                  out <- paste(out, paste0("<i>", i, "</i>"), sep = "<br>")
                }
              }
              
              #otheroutput files
              if(length(other_outputFiles) > 0){
                out <- paste(out, c("<br><b><font color=blue>Other files:</font></b><br>"), sep = "<br><br>")
                
                for(i in other_outputFiles){
                  out <- paste(out, paste0("<i>", i, "</i>"), sep = "<br>")
                }
              }
              
              
              out <- paste0(out, "<br><br>")
              return(out)
            })
            
            clicker <- reactive({
              dat = data.frame(jsonlite::fromJSON(input$clickedData))  
              return(dat)
            })
            
            # if downloads tab is clicked, all output files will get re-written
            observe({
              if(input$inTabset == 'downloads'){
                if(length(input$clickedData) > 0){
                  if(nrow(clicker()) > 0){
                    n.clicker = apply(clicker(), 2, unlist)
                    write.csv(n.clicker, file = paste0(inputFilePrefix, "_", "Clicked_Pathways.csv"))
                  }
                  else{
                    write.csv(data.frame('Clicked'=NA), file = paste0(inputFilePrefix, "_", "Clicked_Pathways.csv"))
                  }
                }
                else{
                  write.csv(data.frame('Clicked'=NA), file = paste0(inputFilePrefix, "_", "Clicked_Pathways.csv"))
                }
                # fwrite(rbindlist(json_1df), file = paste0(inputFilePrefix, "_", "First_Degree_Network.csv"))
                # fwrite(rbindlist(json_2df), file = paste0(inputFilePrefix, "_", "Second_Degree_Network.csv"))
                fwrite(SIGNALoutput.condensed, file = SIGNAL.cond.output.name)
                fwrite(FinalEnrichment.condensed, file = Enrichment.cond.output.name)
                fwrite(Scores_nodes_and_edges, file.name.snae)
                fwrite(non.signal.dat, file = paste0(inputFilePrefix, "_", "nonSIGNALhits.csv"))
                
                ## Update the 'Download' tab
                output$downloadFiles <- renderUI({
                  updateTabsetPanel(session, "inTabset", selected = "downloads")
                  HTML(downloadOut())
                })
              }
            })
            
            
            observeEvent(input$clickedData, {      
              output$ClickedDataTable <- renderDataTable({
                # if(nrow(clicker())==0){
                #   clicker$dat=NULL
                # }
                dat <- datatable(data.frame(clicker()), rownames = TRUE)
                return(dat)
              })
              
            })
            
    
            progress1$inc(1)
            head(E(G))
            
            ###########################################
            ## PathNet / Network Graph Visualization ##
            ###########################################
            
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
      
      ## reactive statement for output files to be printed on download tab
      downloadOut <- reactive({
        outputFiles = list.files(path = './')
        # reordering files so that enrichment and SIGNAL hits files come first and second  and as a seperate list on top
        enrich.spot = grep("SIGNALenrichment", outputFiles)
        hits.spot = grep("_SIGNALhits", outputFiles)
        main_outputFiles = c(outputFiles[c(hits.spot, enrich.spot)])
        
        
        #Select files for non hits
        nonhits.spot = grep("_nonSIGNALhits", outputFiles)
        nonhits_outputFiles = c(outputFiles[c(nonhits.spot)])
        
        #Select files for networks
        network.spot = grep("SIGNALnetwork", outputFiles)
        clickpath.spot = grep("Clicked_Pathways", outputFiles)
        network_outputFiles = c(outputFiles[c(network.spot, clickpath.spot)])
        
        #select all other files not selected above
        
        other_outputFiles = c(outputFiles[-c(enrich.spot, hits.spot, nonhits.spot, network.spot, clickpath.spot)])
        
        
        
        out <- c("<br><b>All files from your SIGNAL analysis for download:</b><br> <br><b><font color=blue>SIGNAL selected hits and pathway enrichments:</font></b><br>")
        
        
        #Main output files
        if(length(main_outputFiles) > 0){
          for(i in main_outputFiles){
            out <- paste(out, paste0("<i>", i, "</i>"), sep = "<br>")
          }
        }
        
        #nonhits output files
        if(length(nonhits_outputFiles) > 0){
          out <- paste(out, c("<br><b><font color=blue>High confidence hits that were not selected by SIGNAL:</font></b><br>"), sep = "<br><br>")
          
          for(i in nonhits_outputFiles){
            out <- paste(out, paste0("<i>", i, "</i>"), sep = "<br>")
          }
        }
        
        #nnetwork output files
        if(length(network_outputFiles) > 0){
          out <- paste(out, c("<br><b><font color=blue>Selected networks and clicked pathway files:</font></b><br>"), sep = "<br><br>")
          
          for(i in network_outputFiles){
            out <- paste(out, paste0("<i>", i, "</i>"), sep = "<br>")
          }
        }
        
        #otheroutput files
        if(length(other_outputFiles) > 0){
          out <- paste(out, c("<br><b><font color=blue>Other files:</font></b><br>"), sep = "<br><br>")
          
          for(i in other_outputFiles){
            out <- paste(out, paste0("<i>", i, "</i>"), sep = "<br>")
          }
        }
        
        
        out <- paste0(out, "<br><br>")
        return(out)
      })
      
      ## Create the 'Download' tab
      output$downloadFiles <- renderUI({
        updateTabsetPanel(session, "inTabset", selected = "downloads")
        
        #outputFiles <- list.files(path = './')
        HTML(downloadOut())
      })

      ## Download all output data
      output$downloadButton <- downloadHandler(
        filename = function(){
          paste("SIGNAL_analysis_output", "zip", sep=".")
        },
        content = function(filename){
          
          outputFiles <- list.files(path = './')
          zip(zipfile=filename, files = outputFiles)
        },
        contentType = "application/zip"
      )
      
      print(paste0("SIGNAL Analysis Time: ", Sys.time() - startA))
      message("Download content completed")
    })
      
      
      ##Download unmapped genes document
      output$unmapped_download <- downloadHandler(
        filename = function() {
          paste("Unmapped_rows,", ".csv", sep="")
        },
        content = function(file) {
          write.csv(data_unmapped, file)
        }
      )
        
        
      #   filename = "unmapped_rows.csv",
      #   content = function(file) {
      #     file.copy(paste0(outDir, "/unmapped_rows.csv"), file)}
      # )
      
      
      ################## User Guide Tab output--------
      output$sampledataset_page <- renderUI({
        tagList(HTML("<p><strong>Guide for using the sample data set for SIGNAL analysis</strong></p>
                     <p>&nbsp;</p>
                     <p><i>Sample dataset:</i></p>"),
                HTML("<br><a href=\"SampleDatasets/LPS_sample_dataset.csv\"><b>Download a Sample SIGNAL Dataset</b></a>
<p>&nbsp;</p>
<p>The sample dataset uses the data from the <em>Sun et al. </em>siRNA study measuring the TNF transcriptional response to LPS in Human Macrophages [1].</p>
<p>&nbsp;</p>
<ol>
<li>To run the analysis in SIGNAL, begin by accessing the SIGNAL interface in your web browser <a href='https://signal.niaid.nih.gov'>https://signal.niaid.nih.gov</a> (Chrome is the recommended browser) and downloading the dataset to your local device.</li>
<li>Click on the &ldquo;Browse&rdquo; icon in the left panel to locate where you have saved the sample file. Select the sample file and click &ldquo;choose&rdquo;.</li>
<li>A dropdown bar with the title &ldquo;Cutoff Type&rdquo; will appear,</li>
</ol>
<p>&nbsp;</p>
<p>There are three ways that a dataset can be split into high-confidence, medium confidence, and low confidence/non hits sets. The sample dataset can be analyzed using either of these approaches.</p>
<ol>
<li>Using a single criterion:
<ol>
<li>Select &ldquo;Zscore&rdquo; from the &ldquo;Cutoff Type&rdquo; dropdown menu.</li>
<li>In the field &ldquo;High Confidence Cutoff Value&rdquo; enter the value -2, in the field below it (&ldquo;Medium Confidence Cutoff Value&rdquo;) enter the value -1.5.</li>
<li>Leave the &ldquo;Add an Additional Criteria&rdquo; unchecked.</li>
</ol>
</li>
</ol>
<p>&nbsp;</p>
<ol start = '2'>
<li>Using two criteria:
<ol>
<li>select &ldquo;Zscore&rdquo; from the &ldquo;Cutoff Type&rdquo; dropdown menu.</li>
<li>In the field &ldquo;High Confidence Cutoff Value&rdquo; enter the value -2, in the field below it (&ldquo;Medium Confidence Cutoff Value&rdquo;) enter the value -1.5.</li>
<li>Click on the &ldquo;Add an Additional Criteria&rdquo; option.</li>
<li>Select &ldquo;OffTarget_pValue&rdquo; from the &ldquo;Column to Use for Secondary Criteria&rdquo; dropdown menu.</li>
<li>Select &ldquo;&ge;&rdquo; from the &ldquo;Direction&rdquo; dropdown menu.</li>
<li>In the field &ldquo;Value&rdquo; enter 0.05.</li>
</ol>
</li>
</ol>
<p>&nbsp;</p>
<ol start = '3'>
<li>Assigning a distinct value to each gene before uploading to SIGNAL and using the assigned value to separate the tiers:
<ol>
<li>select &ldquo;assigned.value&rdquo; from the dropdown menu.</li>
<li>In the field &ldquo;High Confidence Cutoff Value&rdquo; enter the value 1, in the field below it (&ldquo;Medium Confidence Cutoff Value&rdquo;) enter the value 0.5.</li>
<li>Leave the &ldquo;Add an Additional Criteria&rdquo; unchecked.</li>
</ol>
</li>
</ol>
<p>&nbsp;</p>
<p>&nbsp;</p>
<ol start='4'>
<li>Click &ldquo;Analyze my data&rdquo;.</li>
</ol>
<p>&nbsp;</p>
<p>Note: The rest of the guide will follow the analysis using option A (Zscore).</p>
<p>&nbsp;</p>
<p>A progress bar will appear at the bottom right of the browser window. Once it is complete the window will shift to the Enriched Pathways tab listing the enriched pathways identified in the analysis. Gene IDs that were assigned as high confidence hits in the upload file are in blue. Gene IDs that were categorized as medium confidence hits in the input file are in red.</p>
<p>&nbsp;</p>
<p>The names of the enriched pathways can be clicked on to open a new window with a KEGG pathway map overlaying the hits identified in the analysis.</p>
<p>&nbsp;</p>
<ol start='5'>
<li>Click on the &ldquo;Proteasome&rdquo; pathway box. A warning window will popup saying you are leaving an NIH website. Click OK. A new tab will open in your browser with a schematic of the proteasome. The SIGNAL identified hits from the screen are highlighted in blue (high confidence) and in red (medium confidence).</li>
<li>Return to the SIGNAL tab in your browser and click on the &ldquo;Gene Hits&rdquo; tab. The table with the &ldquo;SIGNAL Gene Hits&rdquo; lists all the hits identified by SIGNAL from the uploaded study. Additional columns show which of the other identified hits from the dataset the listed candidate has predicted interactions with and what enriched pathways the interacting genes are members of.</li>
<li>Click on the subtab labeled &ldquo;Graph: Gene Hits by Iteration&rdquo;. A table and graph show the number of high confidence and medium confidence hits selected as candidates through each iteration of SIGNAL. Iteration 0 corresponds to the settings of the upload file with 453 hits described as high confidence. After one iteration a number of the high confidence hits are dropped out and a number of the medium confidence hits are added. The cycle repeats until iteration 3 &amp; 4 that show no changes and thus the iterative analysis terminates.</li>
<li>Click on the &ldquo;High Confidence Hits Not Selected by SIGNAL&rdquo; subtab. All the hits that were designated as high confidence in the upload file but were <em><u>not </u></em>selected as hits by SIGNAL are listed in this table. This helps the user review any of the candidates that were not selected by SIGNAL that the user might want to manually add to the set of selected hits.</li>
<li>Click on the &ldquo;Pathway Enrichments&rdquo; subtab. The enriched pathways with statistics and gene names are shown in a table. The last column includes an &ldquo;Enrich score&rdquo; sorting the enrichment score from highest to lowest, shows a strong enrichment for glycan related pathways.</li>
</ol>
<p>&nbsp;</p>
<p>In addition to the lists of selected hits and pathway enrichments that SIGNAL identifies the platform also enables the exploration of selected hits and possible &lsquo;missing links&rsquo; between enriched pathways. Using a unique and interactive network display, the user can explore predicted interactions between hits associated with specific enrichments and hits outside of the group to infer novel regulatory links.</p>
<p>&nbsp;</p>
<ol start='10'>
<li>Click on the &ldquo;Network&rdquo; tab.</li>
<li>Select (by checking the adjacent box) Pathways 1: Proteasome, 5: Toll-like receptor signaling pathway.</li>
<li>Click &ldquo;&gt;&gt; Create Network Graph&rdquo; at the top of the page. A progress bar will show the progress of the graph. Once complete the &ldquo;Network Path&rdquo; tab will open.</li>
<li>The image can be enlarged by zooming in on the browser image (Command + for Mac, or Ctrl scroll for PC). Adjust the text size using the sliders on the right to make the image more legible and aesthetic.</li>
<li>Hover with your cursor over the names or nodes of different genes to see their interactions and more information in the window at the top right of your browser.</li>
<li>Click on the gene &ldquo;PSMC3&rdquo; that is in the red Proteasome group in the top right of the graph. The predicted interactions will appear, as well as further information in the top right window.</li>
<li>Hover over the highlighted gene &ldquo;TNFRSF1A&rdquo; which is on the left of the graph in the green group. The information window will fill with information about this particular interaction based on the information from the STRING database.</li>
<li>Click on the &ldquo;TNFRSF1A&rdquo; label, the gene ID &ldquo;MAP3K7&rdquo; in the Toll-like receptor signaling pathway group (blue) will highlight among others. Click on &ldquo;MAP3K7&rdquo;.</li>
</ol>
<p>&nbsp;</p>
<p><em>Note: if at any point an unintended or incorrect click was made, the click can be undone by clicking &ldquo;Revert Click&rdquo; at the right side of the browser window. To restart the process, choose &ldquo;Click to reset&rdquo;.</em></p>
<p>&nbsp;</p>
<ol start='18'>
<li>Click &ldquo;Highlight Clicked Pathways&rdquo;. The graph will show the interactive pathway between Proteasome hits and Toll-like receptor signaling pathway hits going through TNFRSF1A.</li>
<li>Click on the &ldquo;Clicked Pathways Table&rdquo; subtab. The pathway clicked through above is presented in table format.</li>
<li>Click on the &ldquo;Download&rdquo; tab. All the analysis generated are available in CSV format and can be download in a zip folder. (Note: the &ldquo;Clicked Pathway Table&rdquo; resets every time the clicked pathway is reset in the graph. To save a specific clicked table before trying a different path click on the download folder to download the files before clicking &ldquo;Click to reset&rdquo; to save the most recent clicked table.) For data security, once the browser window is closed the analysis is deleted from the server.</li>
</ol>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>References</p>
<p>&nbsp;</p>
<ol>
<li>Sun, J., et al., <em>Genome-wide siRNA screen of genes regulating the LPS-induced TNF-alpha response in human macrophages.</em> Sci Data, 2017. <strong>4</strong>: p. 170007.</li>
</ol>
<p>&nbsp;</p>")
        )
      })
      
      
      output$dataPrep_page <- renderUI({
        tagList(HTML(
        "<p><strong>Preparing Your Data:</strong></p>
<p>&nbsp;</p>
<p>There are three requirements for a file to be successfully uploaded for analysis by SIGNAL:</p>
<p>&nbsp;</p>
<ol>
<li>The file must be in a .csv, .txt, or .xlsx format</li>
<li>The file must contain a correctly labeled column with gene IDs (either NCBI EntrezID or HGNC GeneSymbol)</li>
<li>The file must contain a column with numeric values that can be used to separate high confidence hits, medium confidence hits, and non-hits in the data set.</li>
</ol>
<p>&nbsp;</p>
<center><img style='border:3px solid #000000' src='images/SampleInput.png' width='600' ></center>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>A description of these requirements are detailed below:</p>
<p>&nbsp;</p>
<p><strong><u>File format:</u></strong> To upload your data for analysis by SIGNAL first ensure that your document is in .csv, .txt. or .xlsx format.</p>
<p><strong><u>IDs:</u></strong> The file must contain one of the following:</p>
<ol>
<li>A column exactly titled &ldquo;GeneSymbol&rdquo; that has HGNC gene symbols in all the rows</li>
</ol>
<p>or, alternatively,</p>
<ol start='2'>
<li>A column exactly titled &ldquo;EntrezID&rdquo; with NCBI EntrezID in all the rows.</li>
</ol>
<p>Both ID columns, however, do not need to be included in the upload file. If only one is included, SIGNAL will identify the missing ID type and generate a column of the other IDs based on the provided IDs. If both are provided, SIGNAL will use both columns as provided by the input file.</p>
<p><strong><u>Assigning confidence:</u></strong> The upload file must include a column with the numeric values to be used for selecting high and medium confidence hits. The name of this column is up to the user.</p>
<p>&nbsp;</p>
<p>The numeric values can be either continuous values (such as a range of <em>p </em>values or a range of Zscores) or assigned values (such as assigning a value of 1 to all IDs that should be considered &ldquo;high confidence&rdquo;. A value of 0.5 to all IDs that should be considered &ldquo;medium confidence&rdquo; and a value of 0 for all IDs that should be considered &ldquo;non-hits&rdquo;). A more detailed description of how to assign the value and cutoffs to different types of datasets are described (<em>Appendix A: How to segment data into confidence tiers) </em>Each gene ID should only have one value associated with it and be listed in the document only once.</p>
<p>&nbsp;</p>
<p>If the user wishes to use a secondary filtering criteria (such as p-value in addition to fold change), that column must have a different name than the primary column used for assigning confidence cutoffs.</p>
<p>&nbsp;</p>
<p>Additional columns can be included in the document and they will be skipped over by the analysis but preserved in the main download file (<em><u>[name of your upload file]</u></em><em> _SIGNALhits.csv</em>)</p>
        "
        ))
        
      })
      
      output$runAnalysis_page <- renderUI({
        tagList(h4("Running an analysis:", align = 'center'),
                br(),
                br(),
                HTML(
          "<center><img src='images/SideBar_map.png' width='600'></center>
          <p>&nbsp;</p>
<p>Begin by loading the SIGNAL website in your browser using the web address <a href='https://signal.niaid.nih.gov/'>https://signal.niaid.nih.gov/</a></p>
<p>&nbsp;</p>
<p>A panel of dropdown menus will be on the left side of the browser window. Select the parameters that best describe your data and the database setting you want to use for your analysis.</p>
<p>&nbsp;</p>
<center><img style='border:3px solid #000000' src='images/StartAnalysis_page.png' width='600' ></center>
<p>&nbsp;</p>
<p>Some of the parameters come with default options, others require an input from the user. Following is a step-by-step description of the parameter definitions and settings:</p>
<p>&nbsp;</p>
<p><strong><u>Select your organism:</u></strong></p>
<p><strong>&nbsp;</strong></p>
<p>Under this dropdown menu select &ldquo;Human&rdquo; or &ldquo;Mouse&rdquo; based on the gene IDs used in the dataset.</p>
<p>&nbsp;</p>
<p>Default setting: <em>Human.</em></p>
<p>&nbsp;</p>
<p><strong><u>Select a Database for Enrichment Analysis:</u></strong></p>
<p><strong><u>&nbsp;</u></strong></p>
<p>For the enrichment analysis portion, SIGNAL uses the pathway designations curated by the Kyoto Encyclopedia of Genes and Genomes (KEGG). (for information about the KEGG database, see <a href='https://www.kegg.jp/kegg/pathway.html'>https://www.kegg.jp/kegg/pathway.html</a><u>)</u></p>
<p>&nbsp;</p>
<p>Under the <strong>Select a Database for Enrichment Analysis </strong>dropdown, the user can select whether to use only the pathways from the KEGG database that describe biological processes by selecting the <em>KEGG: Biological Processes</em> option, or to use only the pathways associated with disease descriptions by selecting the <em>KEGG: Disease Pathways </em>option. To use both types of pathways the user can select <em>KEGG: All Pathways </em>which includes all of the pathways curated in the KEGG database. &nbsp;</p>
<p>&nbsp;</p>
<p>Default setting: <em>KEGG: Biological Processes.</em></p>
<p>&nbsp;</p>
<p><strong><u>Select Interactions for Network Analysis: </u></strong></p>
<p><strong>&nbsp;</strong></p>
<p>For the network analysis component, SIGNAL uses predicted protein-protein interactions from the STRING database (mapped back to the associated gene name). The default setting is <em>Experimental &amp; Database </em>which incorporates interactions from STRING that have an evidence source in other experimental and curated databases. An additional option is to select <em>Advanced Options </em>from the drop down, once selected a list of six evidence criteria will appear and the user can manually select which interactions to include or exclude based on the evidence criteria of origin. (For a further explanation of each of the possible criteria see <a href='http://string-db.org/cgi/help.pl?sessionId=Z6t6X5gizxo0'>http://string-db.org/cgi/help.pl?sessionId=Z6t6X5gizxo0</a>)<strong>. </strong></p>
<p><strong>&nbsp;</strong></p>
<p>Default setting: <em>Experimental &amp; Database</em></p>
<p>&nbsp;</p>
<p><strong><u>Interaction Confidence for Network Analysis:</u></strong></p>
<p>&nbsp;</p>
<p>The STRING database assigns to each predicted interaction a confidence score ranging from 0 to 1.</p>
<p>&nbsp;</p>
<p>The user can select what confidence score (as defined on the STRING database) to consider for the network interactions used in the analysis. Only interactions whose confidence score cross the selected confidence threshold will be considered in the network analysis portion of the analysis. For simplicity, the provided cutoff options are broken into three groups. <em>Low (&gt; 0.15)</em>, <em>Medium (&gt; 0.4), High (&gt; 0.7). </em>The higher the cutoff, the less predicted interactions are included in the set. The lower the cutoff, the higher the likelihood of false positive interactions being included. Medium confidence is recommended for most cases.</p>
<p>&nbsp;</p>
<p>Default setting: <em>Medium (&gt; 0.4)</em></p>
<p>&nbsp;</p>
<p><strong><u>Choose an input file to upload: </u></strong></p>
<p><strong>&nbsp;</strong></p>
<p>Upload your .csv file by clicking on the &ldquo;Browse&rdquo; button and locating the file in your computer. A progress bar will inform you when the upload is complete. The data from the uploaded file will appear in a table under the &ldquo;Input&rdquo; tab.</p>
<p>&nbsp;</p>
<p><strong><u>Cutoff Type: </u></strong></p>
<p><strong>&nbsp;</strong></p>
<p>Once the file is successfully uploaded a dropdown menu appears under the title &ldquo;Cutoff Type&rdquo;. This is where the numeric column to be used for separating the IDs into different hit confidence groupings will be assigned. The dropdown menu consists of a list of the column names in the uploaded document. Select the column that contains the numeric values to be used for the high confidence/medium confidence cutoffs of your targets.</p>
<p>&nbsp;</p>
<p><strong><u>High Confidence Cutoff Value</u></strong> &amp; <strong><u>Medium Confidence Cutoff Value</u>: </strong></p>
<p><strong>&nbsp;</strong></p>
<p>Type in the numeric values to be used as a cutoff for high confidence and medium confidence hits from your screen. Based on the difference between the high confidence cutoff value entered and Med-confidence cutoff value entered SIGNAL will infer whether the values should be taken as &ldquo;greater than or equal to&rdquo; or &ldquo;less than or equal to&rdquo;.</p>
<p>&nbsp;</p>
<p>(For example if 0.01 is assigned to the <em>High Confidence Cutoff Value </em>and 0.05 is assigned to the <em>Med Confidence Cutoff Value </em>field SIGNAL interprets the input as any ID with a value of 0.01 or less in the selected column should be considered high confidence, any ID with an associated value between 0.01 and 0.05 should be assigned medium confidence, and any ID with an assigned value greater than 0.05 will be assigned as a &ldquo;non-hit&rdquo;. Alternatively, if 2 is assigned to the <em>High-conf Cutoff Value </em>and 1 is assigned to the <em>Med-conf Cutoff Value </em>field SIGNAL interprets the input as any ID with a value of 2 or more in the selected column should be considered high confidence, any ID with an associated value between 1 and 2 should be assigned medium confidence, and any ID with an assigned value less than 1 will be assigned as a &ldquo;non-hit&rdquo;.</p>
<p>&nbsp;</p>
<p><strong><u>Add an Additional Criteria:</u></strong></p>
<p><strong><u>&nbsp;</u></strong></p>
<p>Checking this option allows the user to add an additional condition that all hits must satisfy. This is useful for when an additional metric is to be used, such as <em>p value </em>or <em>cell count</em>, to filter hits in addition to the criteria selected in the &ldquo;Cutoff Type&rdquo; menu.</p>
<p>&nbsp;</p>
<p><strong><u>Column to Use for Secondary Criteria:</u></strong></p>
<p><strong><u>&nbsp;</u></strong></p>
<p>If the &ldquo;Add an Additional Criteria&rdquo; option is selected a drop down menu with all the column names from the uploaded document appears (minus the column name selected under &ldquo;Cutoff Type&rdquo; for the first criteria). Select the column that contains the numeric values that will have to meet the assigned cutoff.</p>
<p>&nbsp;</p>
<p><strong><u>Direction</u></strong> &amp; <strong><u>Value:</u></strong></p>
<p><strong><u>&nbsp;</u></strong></p>
<p>Select the direction (&ge; or &le;) and value for the cutoff to be applied to the secondary criteria. All high-confidence hits and medium confidence hits from the primary criteria that <em>don&rsquo;t</em> meet the cutoff of the secondary criteria are reassigned as &ldquo;non-hits&rdquo;.</p>
<p>&nbsp;</p>
<p><strong><u>Add genome background: </u></strong></p>
<p><strong>&nbsp;</strong></p>
<p>Checking this option adds genes to the list that aren&rsquo;t included in the upload file to be used as a background for statistical enrichment analysis.</p>
<p>&nbsp;</p>
<p>This feature is recommended for when the upload file doesn&rsquo;t include a genome-scale set of &ldquo;no confidence hits&rdquo; (as in when only a partial list of the targets are available and not a list of all targets in the genome). In order to be able to run enrichment analysis statistics on the group of assigned hits it is critical to have a robust background of &ldquo;non-hits&rdquo; against which to measure it. In cases where those IDs are not included, the &ldquo;Add genome background&rdquo; option provides a pseudo genome-scale background for the dataset. If the upload file includes a robust genome-scale background than the &ldquo;Add genome background&rdquo; can remain unclicked and only the IDs in the uploaded file will be used.</p>
<p>&nbsp;</p>
<p>When selected, the added background &ldquo;genes&rdquo; will not appear as suggested hits by the SIGNAL analysis, the background genes are only used as a means to have more robust statistics on the enrichment of pathways. The added background genomes use only the known protein coding genes of the selected organisms (source: Human, HGNC. Mouse, MGI.) that are not in the upload file. SIGNAL uses the difference between the size of selected hits and the number of known protein coding genes for that organism as the number of &ldquo;non-hits&rdquo; for enrichment statistics.</p>
<p>&nbsp;</p>
<p>Default: <em>Unselected.</em></p>
<p><em>&nbsp;</em></p>
<p><strong><u>Optional: Enter gene IDs that should be kept as high confidence hits throughout the analysis:</u></strong></p>
<p><strong><u>&nbsp;</u></strong></p>
<p>This text box can be used to enter any gene IDs (either separated by comma or in separate lines) that should not be filtered out by SIGNAL and should be kept as hits independent of what the iterative analysis finds. This is recommended for when the study identifies genes in a pathway that hasn&rsquo;t yet been annotated by the pathway database or when the user has a particular interest or knowledge about the relevance of these candidates. Adding the gene IDs in this text box will place them in an &ldquo;assigned hit&rdquo; category and their enrichments and network interactions with other hits will be shown in the results.</p>
<p>&nbsp;</p>
<p>It is important here to select to correct gene ID type (EntrezID or GeneSymbol) and to ensure that all the entries are in the correct format.</p>
<p>&nbsp;</p>
<p><strong><u>Analyze my data: </u></strong></p>
<p><strong>&nbsp;</strong></p>
<p>Once the parameters for that analysis have been selected and entered, click the &ldquo;Analyze my data&rdquo; icon and the analysis will begin. A progress bar at the bottom right corner of your browser window will indicate the progress of the analysis.</p>
<p>&nbsp;</p>
<p><strong><u>Reset: </u></strong></p>
<p><strong><u>&nbsp;</u></strong></p>
<p>This tab allows the user to reset all the selections and restart the analysis from the start with default settings. This can be done at any point in the analysis. Clicking the reset tab will remove the uploaded file and all the analysis files generated up to that point. If trying multiple analysis or different settings it is recommended that the analysis be reset between each run by clicking the &ldquo;Reset&rdquo; icon or refreshing the browser.</p>"
        ))
        
      })
      
      output$readResults_page <- renderUI({
        tagList(HTML(
          "<p><strong>Reading SIGNAL Results:</strong></p>
<p><strong>&nbsp;</strong></p>
<p>Once the analysis is complete a range of interactive and downloadable tables are generated to suggest hit selection sets and facilitate exploration of the uploaded data.</p>
<p>&nbsp;</p>
<p>SIGNAL provides analysis output on three levels: selected hits, enriched pathways, and generated networks. The output and data presentation is designed to connect between the three levels of outputs such that, as an example, findings in the selected hits section can be used to prioritize the enriched pathway section, and findings from the enriched pathways data exploration can be used to filter the enriched network, and findings from the generated network can be used to identify critical enrichments.</p>
<p>&nbsp;</p>
<p>The output results are separated into five of the nine visible tabs:</p>
<p>&nbsp;</p>
<p><em>Enriched pathways:</em> An interactive table that lists all the pathway enrichments identified in your data by SIGNAL.</p>
<p>&nbsp;</p>
<p><em>Gene Hits:</em> A number of tables that list the hits selected by SIGNAL and their predicted and associated enrichments.</p>
<p>&nbsp;</p>
<p><em>Network:</em> A table to select up to three pathways from the list of enriched pathways to generate a pathway-network graph.</p>
<p>&nbsp;</p>
<p><em>Network Graph: </em>An interactive figure of the pathway-network graph. These figures will only populate after networks have been selected and plotted from the &lsquo;<em>Network&rsquo;</em> tab.</p>
<p>&nbsp;</p>
<p><em>Download:</em> A list and &ldquo;Download&rdquo; button to locally save all the tables generated in the analysis.</p>
<p>&nbsp;</p>
<p>A description and guide for how to use the information in each tab is below.</p>
<p><strong>&nbsp;</strong></p>
<p><strong><u>Enriched Pathways</u></strong><u>: </u></p>
<p><u>&nbsp;</u></p>
<p>When the analysis is complete, a list of enriched pathways appears in a table under the &ldquo;Enriched Pathways&rdquo; tab. The list includes all pathways that have a <em>p</em> value of 0.05 or less in the completed SIGNAL analysis of the uploaded data. The table also includes a list of the Gene IDs from the pathway that are also hits in the screen. The Gene IDs that were designated as high confidence in the input file are highlighted in blue, Gene IDs that were designated as medium confidence in the input file are highlighted in red.</p>
<p>&nbsp;</p>
<center><img src='images/EnrichPathway_map.png' width='600'></center>
<p>&nbsp;</p>
<p>The table includes columns providing the following information:</p>
<p>&nbsp;</p>
<p><strong>Pathway: </strong>This column includes the names of the enriched pathways. Clicking on the pathway name will open a new tab from the KEGG database showing a schematic of the genes in the pathways with the gene hits from the analysis highlighted. Genes that were marked as high confidence at the start of the analysis are highlighted in blue and those marked as medium confidence are highlighted in red.</p>
<p>&nbsp;</p>
<center><img style='border:3px solid #000000' src='images/KEGGPathway_map.png' width='600' ></center>
<p>&nbsp;</p>
<p><strong>pVal: </strong>The p-values based on a hypergeometric test for the enrichment of each pathway are listed.</p>
<p>&nbsp;</p>
<p><strong>pValFDR: </strong>The p-values with added correction for False Detection Rate are listed.</p>
<p>&nbsp;</p>
<p><strong>pValBonferonni: </strong>The p-values with the Bonferroni corrections for multiple testing are listed.</p>
<p>&nbsp;</p>
<p><strong>TotalGenes: </strong>The total number of genes in the pathway.</p>
<p>&nbsp;</p>
<p><strong>HitGenes: </strong>The number of hit genes as selected by the SIGNAL analysis that are in the pathway.</p>
<p>&nbsp;</p>
<p><strong>HitGeneNames: </strong>The HGNC Gene Symbols of the SIGNAL hit genes in each pathway are listed. Genes that were marked as high confidence at the start of the analysis are highlighted in blue and those that were marked as medium confidence are highlighted in red.</p>
<p>&nbsp;</p>
<p><strong><u>Gene Hits</u></strong></p>
<p><strong><u>&nbsp;</u></strong></p>
<p>The Gene Hits tab is a guide for hit selection from the uploaded data based on the analysis by SIGNAL. The Gene Hits tab is separated into five different tabs</p>
<p>&nbsp;</p>
<p>SIGNAL Gene Hits: <em>A table of all the hits selected by SIGNAL from the high and medium confidence hits in the input file.</em></p>
<p>Gene Hits by Iteration: <em>A table that lists all the hits in the input file along with what category of confidence they were assigned to at each iteration of the SIGNAL analysis.</em></p>
<p>Graph: Gene Hits by Iteration: <em>A table and figure of the number of hits selected at each iteration from the input high and medium confidence categories.</em></p>
<p>High Confidence Hits Not Selected by SIGNAL: <em>A table of hits that were assigned as high confidence in the input file but were not selected as high confidence hits by SIGNAL.</em></p>
<p>Pathway Enrichments: <em>A table of pathway enrichments as assigned by SIGNAL.</em></p>
<p>&nbsp;</p>
<p><strong>SIGNAL Gene Hits: </strong>This table provides a list of hits selected by the SIGNAL analysis. The table also includes supporting information on interacting genes (based on the network criteria that were selected at the start of the analysis) and membership in enriched pathways.</p>
<p>&nbsp;</p>
<center><img src='images/SIGNALHits_map.png' width='600'></center>
<p>&nbsp;</p>
<p>The table&rsquo;s columns and what they include are:</p>
<p>&nbsp;</p>
<p><em>EntrezID: </em>NCBI EntrezID identifier.</p>
<p><em>GeneSymbol:</em> HGNC official Gene Symbol.</p>
<p><em>InputCategory: </em>A designation of &ldquo;HighConf&rdquo; or &ldquo;MedConf&rdquo; based on the original category it was assigned to at the start of the analysis (based on the user provided cutoffs).</p>
<p><em>SIGNALhit: </em>&ldquo;Yes&rdquo; indicates a hit selected by the SIGNAL analysis.</p>
<p><em>Pathway: </em>Names of enriched pathways from the analysis which the associated gene is a member of. (This column only lists pathways that were selected as significantly enriched by the analysis. The individual gene may also be part of other pathways not listed in this table due to the pathway as a whole not being significantly enriched in the final dataset.)</p>
<p><em>InteractingGenes: </em>List of other genes also selected as hits by the SIGNAL analysis that have predicted protein interactions with the assigned gene. The predicted interactions are selected based on the network criteria that was set by the user at the start of the analysis.</p>
<p><em>NetworkGenePathways: </em>List of enriched pathways that the network of interacting genes are individually part of. Number in parenthesis indicates number of interacting genes that are members of each pathway.</p>
<p>&nbsp;</p>
<p><strong>&nbsp;</strong></p>
<p><strong>Gene Hits By Iteration: </strong>This table includes all the columns from the input document with appended columns each iteration of the SIGNAL analysis. The added columns list the confidence category assigned to each gene ID at each analysis step.&nbsp; Columns are designated as either PathwayEnrichment or Network Enrichment followed by the iteration number. Genes counted as high confidence hits in an iteration are indicated as &ldquo;HighConf&rdquo;, genes counted as medium confidence are indicated as &ldquo;MedConf&rdquo;. The &ldquo;SIGNALhit&rdquo; column indicates the gene hits that are counted as final hits in SIGNAL. The top row highlighted in orange indicates the total number of hits considered high confidence at the end of each iteration.</p>
<p>(The table also includes columns with information about the pathways, interactions with other hits, and the pathway membership of the interacting genes as in the SIGNAL Gene Hits table.)</p>
<p><strong>&nbsp;</strong></p>
<p><strong>Graph: Gene Hits By Iteration: </strong>This tab includes a table and graph. The table lists the number of high confidence and medium confidence hits that are selected as hits by SIGNAL at each iteration. (The table starts at iteration 0, corresponding to the input settings. At iteration 0 only the input high confidence hits are considered high confidence, and none of the medium confidence hits are reassigned as high confidence. In subsequent iterations of the analysis some high confidence hits are dropped out and some medium confidence of hits are reassigned as high confidence. The total numbers are listed in the table.</p>
<p>The graph shows the above information in graphical form.</p>
<p>&nbsp;</p>
<center><img src='images/IterationHits_map.png' width='600'></center>
<p>&nbsp;</p>
<p><strong>High Confidence Hits Not Selected By SIGNAL: </strong>This table is a list of all the hits that were not selected as hits by SIGNAL yet were assigned as high confidence hits in the input file. This table can be used to review the hits dropped out by SIGNAL to see if any of them should be manually added to the list of selected hits by SIGNAL based on the user&rsquo;s knowledge and discretion.</p>
<p>&nbsp;</p>
<p><strong>Pathway Enrichments: </strong>This table lists the enriched pathways from the analysis with the associated <em>p</em> values and gene members. This table includes all the information from the Enriched Pathways tab as well as additional columns and details. The pathway enrichment table can be used to further prioritize different subsets of the SIGNAL hit selection and to explore imputed enrichments in the data. An added feature, that can be helpful in the exploration of the data, is the Enrich Score column. The Enrichment Score provides an additional way to prioritize different enrichments beyond the provided <em>p values. </em>The enrichment score represents how strong the enrichment is (that is how many of the genes in the pathway are also in the set of selected hits) and how much of that enrichment is driven by high scoring genes (i.e. of the genes in the hit selection set that drive the enrichment of a particular pathway, how many of them were assigned as &ldquo;high confidence&rdquo; in the input file). The table can be sorted by decreasing Enrichment Scores, the search bar can also be used to look for specific pathway keywords or genes.</p>
<p>&nbsp;</p>
<center><img src='images/PathwayEnrichments_map.png' width='600'></center>
<p>&nbsp;</p>
<p>The columns in the Pathway Enrichments table are:</p>
<p>&nbsp;</p>
<p><em>Pathway: </em>Name of the enriched pathway as it is labeled by the KEGG database.</p>
<p><em>pVal:</em> The p-values based on a two-tailed fisher&rsquo;s exact test for the enrichment of each pathway are listed.</p>
<p><em>pValFDR:</em> The p-values with added correction for False Detection Rate are listed.</p>
<p><em>pValBonferonni:</em> The p-values with the Bonferroni corrections for multiple testing are listed.</p>
<p><em>Genes:</em> The total number of genes in the pathway.</p>
<p><em>HitGenes:</em> The number of hit genes as selected by the SIGNAL analysis that are in the pathway.</p>
<p><em>HighScoreGenes:</em> The number of hit genes as selected by SIGNAL that were also categorized as high confidence based on the user provided cutoff at the start of the analysis</p>
<p><em>HighScoreGeneNames:</em> The HGNC Gene Symbols of the high score hit genes in each pathway.</p>
<p><em>MedScoreGeneNames:</em> The HGNC Gene Symbols of the medium score hit genes in each pathway.</p>
<p><em>EnrichScore:</em> A calculation representing the robustness of the pathway enrichments by the number of genes represented in the SIGNAL dataset and how many of them are high scoring. The EnrichScore is calculated from &nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p><strong><u>Network</u></strong></p>
<p><u>&nbsp;</u></p>
<p>SIGNAL can generate a unique interactive version of an integrated pathway and network graph. The integrated pathway and network figure provides an additional way to explore the enrichments and hits identified by SIGNAL. The user can select up to three pathways from the list of enriched pathways and SIGNAL will generate an interactive graph of all the SIGNAL hits (in the <em>Network Graph</em> tab; see below) that are part of the selected pathways as well as all the SIGNAL hits that have predicted interactions with the &lsquo;pathway&rsquo; genes. This configuration is designed to facilitate an exploration of the enrichment data that goes beyond the curated annotations of the enrichment databases. The interactive and tracking features of the figure can be used to identify possible &ldquo;missing links&rdquo; between different enrichments identified in the dataset. The analysis can also be used to identify novel pathway associations by identifying which targets from the study have predicted interactions with multiple members of an enriched pathway. The various features of this analysis can be used to further prioritize subsets of hits for validation and develop hypotheses for testing.</p>
<p>&nbsp;</p>
<center><img src='images/SelectNetwork_map.png' width='600'></center>
<p>&nbsp;</p>
<p>The names of the pathways are listed in a table. The pathways of interest can be selected by clicking the box beside the name. After the number of pathways have been selected, clicking the &ldquo;Create Network Graph&rdquo; icon at the top will generate the graph and load the next page.</p>
<p>&nbsp;</p>
<p><strong><u>Network Graph</u></strong></p>
<p><strong><u>&nbsp;</u></strong></p>
<p>The network graph page generates an interactive graph of enriched pathways and interacting hits based on the selections from the user in the preceding &ldquo;Network&rdquo; tab. The visual parameters of the graph can be adjusted by the user. The tension of the interaction edges, the text size of the gene IDs, the size of the nodes, and the text of the legend, can all be adjusted using the separate sliders on the right side of the graph.</p>
<p>&nbsp;</p>
<center><img src='images/NetworkGraph_map.png' width='600'></center>
<p>&nbsp;</p>
<p>The design uses specific groupings, color coding, and filtering to make it easier to view and interpret. A more detailed description of these settings is below:</p>
<p>&nbsp;</p>
<p><em>Grouping and colors of hits:</em> The graph separates the groups by color. Gene hits that are members of the selected pathways are shown in blue, red, and brown for the first, second, and third selected pathways, respectively. Gene hits that are annotated as members of more than one pathway selected by the user are placed in separate groups with different colors. Gene hits that are members of pathway 1 and 2 are placed between pathway 1 and pathway 2 genes and highlighted in olive (Hex Color Code #a09d01). Gene hits that are members of pathway 2 and pathway 3 are placed between the pathway 2 and pathway 3 grouped genes and are highlighted in orange (Hex Color Code #ce6702). Gene hits that are members of pathway 1 and pathway 3 are placed to the right of pathway 3 grouped genes and highlighted in saturated dark orange (Hex Color Code #6b5b3e). Gene hits that are annotated as members of all three pathways selected by the user are placed to the left of pathway 1 grouped genes and highlighted in purple (Hex Color Code #6d1c8e). Gene hits that are not annotated as members of any of the three selected pathways, yet have predicted protein-protein interactions with at least one of the hits in one of the pathways, are grouped in the lower half of the graph and highlighted in green. A legend at the top left side of the figure indicates all the relevant color labeling.</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p><em>Filtering:</em> For ease of reading the interactive network graph removes some of the total network and pathway information to make the novel and connecting relationships between hits from different pathways and hits outside the selected pathways easier to highlight. To that end, the network graph does not show predicted interaction between gene hits that are annotated as members of the same pathway. The graph also removes pathway gene hits that don&rsquo;t have any predicted interactions with any of the hits that are outside of the pathway. This configuration is designed to make it easier to explore the interactions suggesting novel connections between enrichment groups and hits from the analysis not annotated within the selected enrichment group.</p>
<p>&nbsp;</p>
<p>The figure includes interactive features that are revealed and tracked based on the user&rsquo;s activity as described below: &nbsp;</p>
<p>&nbsp;</p>
<p><em>Highlighting interactions: </em>Hovering with a cursor over a specific gene (&ldquo;node&rdquo;) highlights all the predicted interactions (&ldquo;edges&rdquo;) of that gene.</p>
<p>&nbsp;</p>
<p><em>Show only highlighted interaction: </em>Clicking on a gene ID or node &lsquo;fixes&rsquo; the highlighted interactions, so that the user can then click on one of the predicted interactions to observe the interactions from the second node.</p>
<p>&nbsp;</p>
<p><em>Evidence source and confidence score of interactions: </em>A panel at the side of the graph provides information about the highlighted interaction, such as the evidence source for the interaction and its confidence score from the STRING database. The panel also lists the number of interactions a node has within the graph. The panel also shows which confidence level the gene ID was categorized in according to the settings at the start of the analysis.</p>
<p>&nbsp;</p>
<p><em>Clicking through the graph to map a novel pathway: </em>In addition to highlighting different interactions, specific interactions can be clicked and resultant pathways mapped.</p>
<p>&nbsp;</p>
<center><img src='images/HighlightNetwork_map.png' width='600'></center>
<p>&nbsp;</p>
<p>After clicking through a string of interacting genes the user can click the &ldquo;Highlight Clicked Pathway&rdquo; icon and all the genes clicked through in the exploration, together with their predicted interactions, are highlighted. &nbsp;The pathways they are members of are also listed in a separate table that can be downloaded with the rest of the analysis at the <em>Download</em> tab.</p>
<p>&nbsp;</p>
<p>To drive further exploration of the data, the SIGNAL platform makes it possible to view which of the gene hits identified by the SIGNAL analysis that are not known members of specific gene sets (&ldquo;Novel&rdquo; genes) have predicted interactions with known members of the gene sets (&ldquo;pathways&rdquo;)</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p>For ease of viewing the network graph is generated with two viewing options and two information tables.</p>
<p>&nbsp;</p>
<p><strong>1<sup>st</sup> Degree Network: </strong>A circular graph showing the hit genes from each selected pathway (separated by group and highlighted by different colors) and only the &ldquo;novel&rdquo; genes from the analysis that have predicted interactions with any of the &ldquo;pathway genes&rdquo;. Hovering with the cursor over a gene name will highlight its interaction path and the name of the gene(s) it is predicted to interact with. As described above, this feature can be used to generate exploratory hypothesis of novel mechanisms and interactions and to identify &ldquo;missing links&rdquo; between predicted biological processes generated by the analysis to be further validated by subsequent research.</p>
<p>&nbsp;</p>
<p><strong>2<sup>nd</sup> Degree Network: </strong>A circular network like the 1<sup>st</sup> Degree Network that also includes &ldquo;novel&rdquo; genes that don&rsquo;t show 1<sup>st</sup> degree connectivity with genes in the pathways of interest but show 2<sup>nd</sup> degree connections to the pathways via predicted interactions with other novel genes identified by SIGNAL that have predicted direct interactions with the pathways of interest. This feature can be used for more complex exploration of possible network reconstruction and to broaden the targets for further exploration.</p>
<p>&nbsp;</p>
<p><strong>Network Table: </strong>&nbsp;The Network table provides the tabulated information for the pathway-based network query selected by the user in the &ldquo;Network&rdquo; panel. The Network table lists all the SIGNAL hit genes that have primary or secondary predicted interactions with SIGNAL hits in the pathways (selected by the user) with the following additional information:</p>
<p>&nbsp;</p>
<center><img src='images/NetworkTable_map.png' width='600'></center>
<p>&nbsp;</p>
<p><em>Group: </em>Which of the (up to) three pathways selected by the user the gene is a member of. If none, the gene is grouped as a &ldquo;Novel&rdquo; interactor with these pathways.</p>
<p><em>Pathway:</em> Which of all of the enriched pathways in the SIGNAL analysis the gene is a member of.</p>
<p><em>Allnet.count: </em>the number of other SIGNAL hit genes the gene has predicted interactions with (based on the user set criteria at the start)</p>
<p><em>Ntwrk.all:</em> The gene names of the genes that have predicted interactions with the individual gene.</p>
<p><em>NtwrkCount.[Name of first selected pathway]: </em>The number of genes in the first selected pathway that are hits by SIGNAL that have predicted interactions with the specific gene target.</p>
<p><em>Ntwrk.[Name of first selected pathway]: </em>The gene names of the genes from the first selected pathways that are hits by SIGNAL that have predicted interactions with the specific gene.</p>
<p><em>NtwrkCount.[Name of second selected pathway]: </em>The number of genes in the second selected pathway that are hits by SIGNAL that have predicted interactions with the specific gene target.</p>
<p><em>Ntwrk.[Name of second selected pathway]: </em>The gene names of the genes from the second selected pathways that are hits by SIGNAL that have predicted interactions with the specific gene.</p>
<p><em>NtwrkCount.[Name of third selected pathway]: </em>The number of genes in the third selected pathway that are hits by SIGNAL that have predicted interactions with the specific gene target.</p>
<p><em>Ntwrk.[Name of third selected pathway]: </em>The gene names of the genes from the third selected pathways that are hits by SIGNAL that have predicted interactions with the specific gene.</p>
<p><em>Total_Path_Hits.net.count: </em>The total number of hits by SIGNAL that are members of all of the selected pathways that have predicted interactions with the specific gene.</p>
<p>&nbsp;</p>
<p><strong>Clicked Pathways Table: </strong>The clicked pathways table highlights the gene (&ldquo;nodes) and interactions (&ldquo;edges&rdquo;) that the user followed and clicked on in the interactive network graph.</p>
<p>&nbsp;</p>
<center><img src='images/SelectedPathway_map.png' width='600'></center>
<p>&nbsp;</p>
<p><em>Name1: </em>Gene Symbol of the clicked node.</p>
<p><em>Node1: </em>Node title (group name followed by Gene Symbols).</p>
<p><em>Parent1: </em>Name of the group the node is in (Pathway name or &ldquo;Novel&rdquo; if it is outside one of the selected pathways.)</p>
<p><em>Interactions1: </em>Number of predicted interactions between the node and other hits selected by the analysis.</p>
<p><em>Screen.Input: </em>The confidence category this gene was assigned to at the start of the analysis.</p>
<p><em>Name2, Node2, Parent2, Interactions2: </em>Gene symbol, node title, group, and number of predicted interactions of the second node.</p>
<p><em>Weight: </em>Assigned confidence weight for the predicted interaction in the STRING database.</p>
<p><em>Source: </em>Evidence source for the predicted interaction in the STRING database.</p>
<p><em></em></p>
<p>&nbsp;</p>
<p>Both graphs can be viewed as interactive HTMLs and screen grabbed for saving the image. .csv files of the Network Table and Clicked Pathway Table can be downloaded from the &ldquo;Download&rdquo; tab.</p>
<p>&nbsp;</p>
<p><em>Note: the &ldquo;Clicked Pathway Table&rdquo; resets every time the clicked pathway is reset in the graph. To save a specific clicked table before trying a different path click on the download folder to download the files before clicking &ldquo;Click to reset&rdquo; to save the most recent clicked table. </em></p>"
        ))
      })
      
      output$saveAnalysis_page <- renderUI(
        tagList(
          HTML(
            "<p><strong>Saving and securing your analysis:</strong></p>
<p>&nbsp;</p>
<p>The Download tab contains a list of all the files generated by the analysis and available to download. The &ldquo;Download all files&rdquo; icon downloads a zipped folder of all the listed files in .csv format.</p>
<p>&nbsp;</p>
<p>As more analysis are added within the same session, the analysis files are added to the zipped folder for download. Once the session ends the analysis files are deleted.</p>
<p>&nbsp;</p>
<p>All the analysis download files begin with the same file name as the file uploaded by the user at the beginning of the analysis with the different titles appended to the end.</p>
<p>&nbsp;</p>
<center><img style='border:3px solid #000000' src='images/Download_map.png' width='600' ></center>
<p>&nbsp;</p>
<p>The file names and their contents are:</p>
<p>&nbsp;</p>
<p><strong>[name of input file]_SIGNALhits.csv: </strong>A table listing all the hits by SIGNAL analysis with enriched pathway and interacting network information.</p>
<p>&nbsp;</p>
<p><strong>[name of input file]_SIGNALenrichment.csv: </strong>A table of enriched pathways identified by the analysis with their associate p-values, FDR values, Bonferroni values, and high confidence and medium confidence hit genes included.</p>
<p>&nbsp;</p>
<p><strong>[name of input file]_ nonSIGNALhits.csv: </strong>a table listing all the genes that were assigned as high confidence hits at the start of the analysis (by the user provided criteria) but were not selected as hits by the end of the SIGNAL analysis.</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p><strong>[name of input file]_SIGNALnetwork_[Name(s) of selected pathways].csv:</strong> A table of all the hit genes selected by SIGNAL analysis with their predicted membership or interactions with genes in the pathways selected under the &ldquo;Network&rdquo; tab.</p>
<p>&nbsp;</p>
<p><strong>[name of input file] _Clicked_Pathways.csv:</strong> A table of all the genes clicked on by the user in the &ldquo;Highlight Clicked Network&rdquo; feature under the network tab.</p>
<p>&nbsp;</p>
<p><strong>Unmapped_rows: </strong>A table of all the rows with GeneSymbols that couldn&rsquo;t be mapped to EntrezIDs. (This document only shows up when the upload file only included a GeneSymbol column and EntrezID mapping was done by SIGNAL).</p>"
          )
        )
      )
      
      output$segmentData_page <- renderUI(
        tagList(
          HTML("<p><strong>Appendix A: How to segment data into confidence tiers</strong></p>
<p><strong>&nbsp;</strong></p>
<p>There are three ways that a dataset can be split into high-confidence, medium confidence, and low confidence/non hits sets. Different datasets require different approaches and different ratios.</p>
<p>&nbsp;</p>
<p>While there&rsquo;s no rule as to what a good number of candidates is to be assigned to each of the data tiers, in our tests we have found that a 1:2 ratio of high confidence hits to medium confidence hits (with all other candidates serving as background/non-hits) work best (i.e. ~400 high confidence hits, ~800 medium confidence hits).</p>
<p>&nbsp;</p>
<p>(The pathway enrichment results from the analysis can serve as &lsquo;gut check&rsquo; for whether some of the expected enrichments show up. The rapid speed with which each analysis is completed on the platform makes it possible to try different cutoffs and different approaches to see which gives a recognizable yet informative result.)</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<p><strong>Data segmentation can be achieved using a single readout, dual measurements, or multiple assays.&nbsp; </strong></p>
<p>&nbsp;</p>
<p>Data can be segmented into confidence scores either by using the normalized readout from an assay or by combining readouts and studies.</p>
<p>&nbsp;</p>
<p>Below is a guide for how to set cutoffs and assign confidence to tiers using three different approaches:</p>
<p>&nbsp;</p>
<p>&nbsp;</p>
<ol style='list-style-type: upper-alpha;'>
<li><strong>Using a single readout or value:</strong></li>
</ol>
<p>&nbsp;</p>
<p><em>This approach is most recommended when using readouts like Z score or Expression that have already been normalized and corrected for outliers.</em></p>
<p><em>&nbsp;</em></p>
<center><img src='images/Segment_1criteria.png' width='600'></center>
<p><em>&nbsp;</em></p>
<p>To segment data into three tiers based on a single readout, simply assign those values to a specific column in your input file. Choose a &ldquo;stringent&rdquo; cutoff for the &ldquo;High Confidence Cutoff Value&rdquo; in SIGNAL and a more &ldquo;lenient&rdquo; cutoff for the &ldquo;Medium Confidence Cutoff Value&rdquo; in SIGNAL.</p>
<p>&nbsp;</p>
<ol style='list-style-type: upper-alpha;' start = '2'>
<li><strong>Using two readouts or values:</strong></li>
</ol>
<p>&nbsp;</p>
<p><em>This approach is most recommended when using values like Fold Change where an additional readout such as the replicate p Value is helpful to increase confidence in the selected genes.</em></p>
<p><em>&nbsp;</em></p>
<center><img src='images/Segment_2criteria.png' width='600'></center>
<p><em>&nbsp;</em></p>
<p>To segment data into three tiers based on two readouts, choose the readout that will be broken into a &ldquo;stringent&rdquo; cutoff and a &ldquo;lenient&rdquo; cutoff (Fold Change for example) assign those values to a specific column in your input file. Choose a &ldquo;stringent&rdquo; cutoff for the &ldquo;High Confidence Cutoff Value&rdquo; in SIGNAL and a more &ldquo;lenient&rdquo; cutoff for &ldquo;Medium Confidence Cutoff Value&rdquo; in SIGNAL.</p>
<p>Assign the second set of values to a different column in your input file. Choose a cutoff which all high confidence and medium confidence hits must meet and enter those values into the fields under &ldquo;Add an Additional Criteria&rdquo;.</p>
<p>&nbsp;</p>
<ol style='list-style-type: upper-alpha;' start = '3'>
<li><strong>Assigning a distinct value to each gene before uploading to SIGNAL and using the assigned value to separate the tiers:</strong></li>
</ol>
<p><em>&nbsp;</em></p>
<p><em>This approach is most recommended when using multiple criteria or when combining the results of multiple studies. While the SIGNAL sidebar menu doesn&rsquo;t have a direct option to combine these studies the information can be combined outside of SIGNAL (such as in Excel or in R) and then run as a SIGNAL analysis.</em></p>
<p><em>&nbsp;</em></p>
<center><img src='images/Segment_3criteria.png' width='600'></center>
<p><em>&nbsp;</em></p>
<p>In Excel or in R, create a new column in the file you upload. In the new column, assign a value of 1 to gene candidates that meet the high confidence criteria (i.e. the gene candidate is a hit in more than one of the studies or in all the studies being compared). Then assign a value of 0.5 in the same column to all the gene candidates that meet the criteria for medium confidence (i.e. the gene candidate is a hit in just one or two of the studies being compared.). Assign a value of 0 in the column to all the other gene candidates that don&rsquo;t meet any of the criteria you&rsquo;ve set.</p>
<p>&nbsp;</p>
<p>A similar approach can be used when selecting hits based on multiple guides in a gene perturbation screen.</p>
<p>&nbsp;</p>
<p>Once the values have been added upload your file to SIGNAL. Under &ldquo;Cutoff Type&rdquo; select the name of the column in which you have added the new values. In the field &ldquo;High Confidence Cutoff Value&rdquo; enter the value 1, in the field below it (&ldquo;Medium Confidence Cutoff Value&rdquo;) enter the value 0.5. This will run the iterative analysis using the high, medium, and low confidence tiers that you&rsquo;ve assigned outside of SIGNAL.</p>"
            
          )
        )
      )
      
      
      output$dataSecurity_page <- renderUI(
        tagList(
          HTML("<p><strong>Appendix B: Data security</strong></p>
<p>The SIGNAL application web interface is hosted by the <a href= 'https://www.niaid.nih.gov/about/cyber-infrastructure-computational-biology-contacts'>National Institute of Allergy and Infectious Disease (NIAID) Office of Cyber Infrastructure and Computational Biology (OCICB).</a>. The analysis of uploaded data is run behind two internet security firewalls. Access by external users to the SIGNAL interface pass through two secure firewalls. The incoming request first passes through the NIH web hosting firewall after which the analysis and SIGNAL application go through the NIAID firewall where the analysis is hosted. Migrations to external websites outside of the firewall (i.e. the KEGG interface) are accompanied with a warning message for the user.</p>
<p>SIGNAL uses a secure encrypted HTTPS connection, and all requests are handled using encrypted connections. Using these encrypted connections, only the browser from the IP address where the request originated from can access the data generated and uploaded. When a file is uploaded, a new unique directory is created where the input file and all the subsequently generated analysis files are temporarily saved, ensuring that only the user generating the connection can access the directory with their files. The directory is kept on the server only for the duration of the session (i.e. as long as the user is using the site). Once the sessions ends (i.e. close of browser window or move to a new site) the directory and all its files are removed from the SIGNAL server. This decreases the security risk for the user and ensures that the results are only stored locally after the analysis.</p>
<p>Data collected for each session is limited to the country where the request comes from and the time spent using the site (and specific pages). File names, analysis choices, user IDs, and results are neither collected nor stored.</p>
<p>An additional option to increase the security of the data and generated analysis is to download SIGNAL source code and run the application locally on your own device see 'Appendix C: Running SIGNAL with alternative or bespoke databases'.</p>
<p>For any further questions about data security, please reach out to us via the &lsquo;Contact us&rsquo; tab under the &ldquo;Help&rdquo; tab.</p>"
            
          )
        )
      )
      
      
      output$standaloneR <- renderUI(
        tagList(
          HTML("<p><strong>Appendix C: Running SIGNAL with alternative or bespoke databases</strong></p> <p>&nbsp;</p>"),
          HTML( "<p>To make SIGNAL an adaptable framework for iterative analysis with different datasets and databases beyond the databases and settings used on this platform, am R script version of a standalone SIGNAL function can be downloaded. The SIGNAL function relies on calling two separate analysis function, a pathway enrichment function and a network analysis function. The master SIGNAL function applies the pathway and network function iteratively, and the results are tested for when the analysis converges on a single set.</p>"),
          HTML("<br><a href=\"RscriptDownload/SIGNAL_R_function.R\"><b>Download SIGNAL as a standalone R function</b></a><br><br><br>
<p>The list of input variables that can be selectively assigned in the adaptable SIGNAL function in R and their required formats are:</p>
<p>&nbsp;</p>
<p><em>screen.datafame: </em>A data frame of the screen.</p>
<p><em>ID.column:</em> A column within the screen.dataframe for the identifiers of the targets (EntrezID, GeneSymbol, etc.).</p>
<p><em>criteria.column:</em> A column within the screen.dataframe of the criteria for being considered a hit.</p>
<p><em>highconf.criteria:</em> A criteria each target has to meet to be considered a 'high confidence' hit.</p>
<p><em>midconf.criteria: </em>A criteria each target has to meet to be considered a 'mid confidence' hit.</p>
<p><em>criteria.setting:</em> Whether the function should be using 'equal', 'greater than or equal', or 'less than or equal' when assessing if confidence criteria are met. criteria.setting input should be in the format of 'equal', 'greater', or 'less'.</p>
<p><em>enrichment.dataframe:</em> A data frame to be used for pathway membership in the format of a column of IDs (should be same as ID column in screen.dataframe in ID type and column title) and a column of which group they are part of (each ID~group relationship needs to be in its own separate row).</p>
<p><em>enrichment.title:</em> Name of the column with the names of the enrichment groups the targets are members of.</p>
<p><em>stat.test:</em> Name of the statistical test to be used for measuring enrichment confidence. Needs to be in the format of either 'pVal', 'FDR', or 'Bonferroni'.</p>
<p><em>test.cutoff:</em> A numeric value which a less than value in stat.test will be considered a significant enrichment.</p>
<p><em>network.igraph:</em> an igraph of the network to be used for network analysis (network igraph must use the same ID type as screen.dataframe)</p>
<p>&nbsp;</p>
<p>The user provided variables are then used to apply the iterative function as in the previous paragraph. The adaptable version of SIGNAL broadens the scope of its application beyond the use of the specific databases and settings it was designed in.</p>
<p>&nbsp;</p>
<p>The SIGNAL function provides an output in the format of a R script list that contains three data frames:</p>
<ol>
<li>The input data frame with an appended 'SIGNAL.hit' column.</li>
<li>A data frame of high confidence and medium confidence designation at each iteration of the analysis.</li>
<li>A data frame of final SIGNAL enrichments from the provided enrichment data frame.</li>
</ol>"
        ))
      )
      
      ####Downlaod button for SIGNAL R script
      output$sampleDataset <- downloadHandler(
        filename = "SIGNAL_R_function.R",
        content = function(file) {
          file.copy("www/RscriptDownload/SIGNAL_R_function.R", file)}
      )

      
      
    # Contact us by sending email to signal-team@nih.gov
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
        # Send the email to SIGNAL team
        send.mail(from = input$from,
                  to <- c("signal-team@nih.gov"),
                  subject = input$subject,
                  body = paste(paste0("This email is from: ", input$userName, " [", input$from, "]"), input$message, sep="<br/>"),
                  smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = "signal.lisb@gmail.com", passwd = "LISB@NIH", ssl = TRUE),
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
        src="SIGNAL_userguide_V4.pdf",
        height=700,
        width=800
      )
    })
    
    ## Change log
    output$changeLog <- renderUI({
      tagList(
        HTML("Change Log - version and update history"),
        br(),
        br(),
        HTML(paste0(SIGNAL.version)),
        br(),
        HTML(paste0("Latest update: ", SIGNAL.update.date, ".")),
        br(),
        br(),
        br(),
        HTML("Database Updates Log:"),
        br(),
        HTML(paste0("KEGG Pathways: version ", KEGG.version, ". Latest update: ", KEGG.update.date, ".")),
        br(),
        HTML(paste0("STRING Network: version ", STRING.version, ". Latest update: ", STRING.update.date, "."))
        )
      
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

