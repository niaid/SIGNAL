library(shiny)
library(DT)

server <- function(input, output, session) {
  X = data.frame(
    ID = c("Click here then Jump to tab2 and pass x=1 and y=10 to tab2", 
           "Click here then Jump to tab2 and pass x=2 and y=2 to tab2",
           "Click here then Jump to tab2 and pass x=3 and y=4 to tab2"),
    x = c(1,2,3),
    y = c(10,2,4)    
  )
  
  output$datatable = renderDataTable({X}, selection = "single", server = FALSE,
                                     options = list(paging=FALSE,
                                                    searching=FALSE,
                                                    filtering=FALSE,
                                                    ordering=FALSE)
  )
  observeEvent(input$datatable_rows_selected, {
    row <- input$datatable_rows_selected 
    output$text <- renderText({paste("X =", X[row, "x"], "Y =", X[row, "y"])})
    updateTabsetPanel(session, "mainPanel", selected = "tab2")
  })  
}

ui <- fluidPage(
  
  tabsetPanel(id = "mainPanel", 
              tabPanel("tab1",dataTableOutput("datatable")),
              tabPanel("tab2",textOutput("text"))
  )
)

shinyApp(ui = ui, server = server)