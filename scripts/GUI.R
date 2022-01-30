library(shiny)
source("main.R")

# User Interface
ui <- fluidPage(
  
  titlePanel("Conditions"),
  
  sidebarLayout(
    position = "left",
    sidebarPanel(
      selectInput("file",
                  "View comparison",
                  choices = c(
                    "Lactose vs Glucose 24h",
                    "Lactose vs Glucose 48h",
                    "90% Glucose, 10% Lactose vs Glucose 24h",
                    "90% Glucose, 10% Lactose vs Glucose 48h",
                    "75% Glucose, 25% Lactose vs Glucose 24h",
                    "75% Glucose, 25% Lactose vs Glucose 48h"
                  )
      ),
      checkboxInput("all",
                    "View all comparisons"
      ),
      radioButtons("radio", h3("Radio buttons"),
                   choices = list("Choice 1" = 1, "Choice 2" = 2,
                                  "Choice 3" = 3),selected = 1),
    ),
    
    mainPanel(
      plotOutput("distPlot"),
      #plotOutput("allPlot")
    )
  )
)

# Server
server <- function(input, output) {
  
  output$distPlot <- renderPlot({
    if(!input$all){
      sub_table <- up_down_table[which(up_down_table$label == input$file),]
      plot_diff_results(sub_table, fontsize = 12)
    }
    else{
      plot_diff_results(up_down_table, fontsize = 6)
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
