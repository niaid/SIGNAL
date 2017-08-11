#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(readr)
source("./Rscripts/*.R")

# Define UI for application that draws a histogram
ui <- fluidPage( 
  
   # style
   theme = "./css/triage.css",

   # Application title
   headerPanel(includeHTML("header.html")),

   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
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
         fileInput(inputId='file1', 'Choose an inputfile to upload'
                     ),
         submitButton("Analyze my data!", icon("refresh")),

         width = 3
      #   tags$hr(),
      #   checkboxInput('header', 'Header', TRUE),
      #   radioButtons('sep', 'Separator',
      #                c(Comma=',',
      #                  Semicolon=';',
      #                  Tab='\t'),
      #                ','),
      #   radioButtons('quote', 'Quote',
      #                c(None='',
      #                  'Double Quote'='"',
      #                  'Single Quote'="'"),
      #                '"'),
      #   tags$hr(),
      #   p('If you want a sample .csv or .tsv file to upload,',
      #     'you can first download the sample',
      #     a(href = 'mtcars.csv', 'mtcars.csv'), 'or',
      #     a(href = 'pressure.tsv', 'pressure.tsv'),
      #     'files, and then try uploading them.'
      #   )
      ),
      # Show a plot of the generated distribution
      mainPanel(
         textOutput("title"),
         tableOutput("contents")
      )
   ),
   # Show a footer using the header style
   headerPanel(includeHTML("footer.html"))
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$contents <- renderTable({
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    data <- read.csv(inFile$datapath)
    
    # display the 
    inputDataSummary <- c("Your input data have ",  nrow(data),  "rows and ",  ncol(data), "columns!")
    output$title <- renderText(inputDataSummary)
    
    # display the input data in a table
    data
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

