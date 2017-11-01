## app.R ##
library(shiny)
library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(title = "<b>TRIAGE - Throughput Ranking by Iterative Analysis of Genomic Enrichment", 
                  titleWidth = 400),
  dashboardSidebar(),
  dashboardBody()
)

server <- function(input, output) { }

shinyApp(ui, server)