#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#library(shiny)
#library(shinydashboard)
#library(readr)

if (interactive()) {
# Define UI for application that draws a histogram
ui <- fluidPage( 
  # style
  theme = "./css/triage.css",
  
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
      selectInput("cutoff_type", "Cutoff type", c("P value", "FDR", "Z score")
      ),
      # cutoff values depending the cutoff method chosen
      uiOutput("cutoffValues"),
      actionButton("goButton", "Analyze my data !", icon("angle-double-right"), 
                   style="padding:4px; font-size:120%; color: #fff; background-color: rgb(1, 81, 154); border-color: #2e6da4"),
      width = 3
    ),
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(id = "inTabset",
        tabPanel(title = "Input", value = "contents", 
                 textOutput("dataSummary"), 
                 textOutput("dataHead"), 
                 tableOutput("contents"), 
                 textOutput("organism"),
                 textOutput("pathway"),
                 textOutput("network"),
                 textOutput("cutoff_type"),
                 textOutput("cutoff_value")
        ),
        tabPanel(title = "Status", value = "status", verbatimTextOutput("status")),
        tabPanel(title = "Result", value = "results", tableOutput("results"))
      )    
    )
  ),

  # Show a footer using the header style
  headerPanel(includeHTML("footer.html"))
)


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  output$cutoffValues <- renderUI({
    if (is.null(input$cutoff_type))
      return()

    switch(input$cutoff_type,
           "Z score" = selectInput("cutoff_value", "Cutoff value",
                                   choices = c("+/-1.96", "+/-2.58"),
                                   selected = c("+/-1.96")
           ),
           "P value" = selectInput("cutoff_value", "Cutoff value",
                                   choices = c("<0.05", "<0.01"),
                                   selected = "<0.05"
           ),
           "FDR" = selectInput("cutoff_value", "Cutoff value",
                               choices = c("<5%", "<1%"),
                               selected = c("<5%")
           )
      )
  })

  # Read in the input fie  
  output$contents <- renderTable({
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    data <- read.csv(inFile$datapath)
    
    # display the input file dimension
    output$dataSummary <- renderText(c("Your input data have ",  nrow(data),  "rows and ",  ncol(data), "columns."))

    # display the input data in a table
    output$dataHead <- renderText("Here are the first 10 rows:")
    head(data, 10)
  })
  
  # parameters
  output$organism <- eventReactive(input$goButton, paste("Organism:", input$organism))
  output$pathway  <- eventReactive(input$goButton, paste("Pathway:", input$pathway))
  output$network  <- eventReactive(input$goButton, paste("Network:", input$network))
  output$cutoff_type <- eventReactive(input$goButton, paste("Cutoff Type:", input$cutoff_type))
  output$cutoff_value <- eventReactive(input$goButton, paste("Cutoff Value:", input$cutoff_value))
  
  # Upon job submission
  output$status <- eventReactive(input$goButton, 'Job submitted!')
  
  observeEvent(input$goButton, {
    updateTabsetPanel(session, "inTabset",
                      selected = "status")
  
      #source("~/TRIAGE/Rscripts/PathNet-iteration_V3.R")
      source("~/TRIAGE/Rscripts/pathway_iteration.R")
      source("~/TRIAGE/Rscripts/Network_iteration_V3.R")
      
      # Global variables
      scriptDir <- "~/TRIAGE/Rscripts/"
      inputDir <- "~/TRIAGE/inputOutputs/TRIAGEinputFiles/"
      outputDir <- "~/TRIAGE/inputOutputs/TRIAGEoutputFiles/"
      dataDir <- "~/data/"
      source(paste0(scriptDir, "pathway_iteration.R"))
      Organism <- input$organism
      network.type <- "hSTRINGppi.hi"
      
      outDir <- "~/TRIAGE/inputOutputs/TRIAGEoutputFiles"
      outputFileName <- paste0("TRIAGEoutput_HuTNF_CSAfdr_top5percCUTOFF_KEGG_", network.type, ".csv")

      # 1) Seed Pathway Analysis
      setwd(dataDir)
      #pathway.types <- c("KEGG", "Reactome", "Gene_Ontology")
      #pathway.type <- pathway.types[1]
      pathway.type <- input$pathway
      
      pathwayData <- read.csv(file = paste0(dataDir, "Pathways/", pathway.type, Organism, ".csv"))
      
      #seedName <- "IAMinput_HuTNF_manual.csv"
      #seedName <- "TRIAGEinput_HuTNF_CSAfdr_5percCO.csv"
      #setwd(inputDir)
      
      #siRNA.Score <- read.csv(seedName, stringsAsFactors = F)
      siRNA.Score <- read.csv((input$file1)$datapath, stringsAsFactors = F)
      proxyScore <- "Replicate1" 
      iteration <- 1
      counter <- TRUE
      
      while (counter == TRUE) {
        
        Hits <- siRNA.Score$EntrezID[siRNA.Score[[proxyScore]] == 1]
        nonHits <- setdiff(siRNA.Score$EntrezID, Hits)
        
        outPrefix <- paste(pathway.type, iteration, sep = "_")  
        
        # 1) Contraction - [Pathway Analysis]
        siRNA.Score <- ComputeEnrichment(pathwayData, Hits, nonHits, outPrefix, siRNA.Score, iteration)
        siRNA.Score <- data.frame(siRNA.Score, temp = 0, stringsAsFactors = FALSE)
        kName1 <- paste0("KEGG.class.iteration", iteration)
        kName2 <- paste0("KEGG.", iteration)
        names(siRNA.Score)[names(siRNA.Score) == "temp"] <- kName1
        siRNA.Score[[kName1]][siRNA.Score$KEGG == "Yes" & siRNA.Score[[proxyScore]] > 0] <- 1
        siRNA.Score[[kName1]][siRNA.Score$KEGG != "Yes" & siRNA.Score[[proxyScore]] > 0] <- 0.5
        names(siRNA.Score)[names(siRNA.Score) == "KEGG"] <- kName2
        
        hit.Genes <- siRNA.Score$EntrezID[siRNA.Score[[kName1]] == 1]
        
        # 2) Expansion - [Network Analysis]
        system.time(
          source(paste0(scriptDir, "Network_iteration_V2.R"))
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

        #print(paste("iteration: ", iteration))
        cat(paste("iteration: ", iteration, "\n"), file=stderr())
        
        proxyScore <- nName2
        iteration <- iteration + 1
      }
      
      ### Append Enrichment Info --------------------------------------------
      pval_threshold <- 0.05
      iterationNum = iteration - 1
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
      
      ## Save the result into an output file in the TRIAGEoutputFiles folder
      out <- merge(siRNA.Score, pathDF, by = "GeneSymbol", all = T)
      setwd(outDir)
      write.csv(out, file = outputFileName, row.names = F)
      completed <- TRUE
      
      ## Display partial result in the 'results' tab
      observe({
        if(completed) {
          updateTabsetPanel(session, "inTabset",
                            selected = "results")
          
          # Display the complete table
          #output$results <- renderTable({out[100:115,]})
          # Display the pathway table
          output$results <- renderTable({head(pathEnrich)})
        }
      })
    })
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server)
}