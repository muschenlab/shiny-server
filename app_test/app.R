library(shiny)
library(reactlog)

options(shiny.reactlog=TRUE)

dt <- mtcars

ui <- fluidPage(
  # sliderInput(inputId = "num",
  #             label = "choose a number",
  #             value= 20, min = 1 , max = 100),
  selectInput(inputId = "lst",
              label = "select which variable to plot",
              choices = colnames(dt)),
    

  plotOutput(outputId = "hist")
  
)

server <- function(input,output,session){
  
  
  output$hist <- renderPlot({hist(dt[,input$lst],main = paste0(input$lst) )})
  
}


shinyApp(ui = ui, server = server)