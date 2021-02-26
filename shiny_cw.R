# shiny app for interactive inspection of data

library(shiny)
source("cw_functions.R")

ui <- fluidPage(
    
    titlePanel("Shiny Cognition Wall"),
    
    sidebarLayout(
        sidebarPanel(
            
        ),
        mainPanel(
            
        )
    )
)

server <- function(input, output) {
    
}

shinyApp(ui, server)
