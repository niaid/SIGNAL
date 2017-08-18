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
#source("~/TRIAGE/Rscripts/PathNet-iteration_V2.R")
#source("~/TRIAGE/Rscripts/Network_iteration_V2.R")
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
                  choices = c("STRING", "BioGRID") 
      ),
      fileInput(inputId= "file1", 
                label = 'Choose an inputfile to upload'
      ),
      selectInput("cutoff_type", "Cutoff type", c("P value", "FDR", "Z score")
      ),
      # cutoff values depending the cutoff method chosen
      uiOutput("cutoffValues"),
      actionButton("goButton", "Analyze my data !", icon("angle-double-right"), 
                   style="padding:4px; font-size:120%; color: #fff; background-color: rgb(1, 81, 154); border-color: #2e6da4")
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
        tabPanel(title = "Log", value = "logs", verbatimTextOutput("logs")),
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
           "P value" = selectInput("cutoff_value", "Cutoff value",
                                   choices = c("<0.05", "<0.01"),
                                   selected = "<0.05"
           ),
           "FDR" = selectInput("cutoff_value", "Cutoff value",
                               choices = c("<5%", "<1%"),
                               selected = c("<5%")
           ),
           "Z score" = selectInput("cutoff_value", "Cutoff value",
                                   choices = c("+/-1.96", "+/-2.58"),
                                   selected = c("+/-1.96")
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
    output$dataSummary <- renderText(c("Your input data have ",  nrow(data),  "rows and ",  ncol(data), "columns!"))

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
  output$logs <- eventReactive(input$goButton, 'Job submitted!')
  
  observeEvent(input$goButton, {
    updateTabsetPanel(session, "inTabset",
                      selected = "logs")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
}