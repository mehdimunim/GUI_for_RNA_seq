library(shiny)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("main.R")

# User Interface
ui <- fluidPage(tabsetPanel(
  tabPanel("View RPKMs",
           sidebarLayout(
             sidebarPanel(
               selectInput("analysis.type",
                           "View analysis",
                           choices = analysis$analysis)
             ),
             mainPanel(plotOutput("plotRPKM"))
           )),
  
  
  tabPanel(
    "Condition comparison",
    sidebarLayout(
      sidebarPanel(
        selectInput("comparison",
                    "View comparison",
                    choices = analysis$analysis),
        checkboxInput("allComparison",
                      "View all comparisons")
      ),
      mainPanel(plotOutput("plotComparison"))
    )
  ),
  
  tabPanel(
    "Clustering",
           sidebarLayout(
             sidebarPanel(
               selectInput(
               "resultType",
               "Result type",
               choices = c(
                 "Number of genes" = "num",
                 "Expression profile" = "expProf",
                 "Average expression profile" = "avgExp"))
             ),
             mainPanel(plotOutput("plotClustering"))
           )),
  
  tabPanel("Test",
           sidebarLayout(
             sidebarPanel(
               selectInput("analysisTest",
                           "View analysis",
                           choices = analysis$analysis)
             ),
             mainPanel(verbatimTextOutput("printTest"))
           )),
  
  tabPanel("About",
           fluidRow(column(
             6,
             includeMarkdown("../README.md")
           )))
  
))

# Server
server <- function(input, output) {
  output$plotRPKM <- renderPlot({
    x <- analysis[which(analysis$analysis == input$analysis.type), ]
    plot_lfpkm_comparison(unlist(x$fpkm1), unlist(x$fpkm2), x$analysis)
  })
  
  output$plotComparison <- renderPlot({
    if (!input$allComparison) {
      sub_table <-
        up_down_table[which(up_down_table$label == input$comparison), ]
      plot_diff_results(sub_table, fontsize = 12)
    }
    else{
      plot_diff_results(up_down_table, fontsize = 6)
    }
  })
  
  output$printTest <- renderPrint({
    x <- analysis[which(analysis$analysis == input$analysisTest), ]
    summarize_test(x$fpkm1, x$fpkm2, x$analysis, alpha)
  })
  
  output$plotClustering <- renderPlot({
    if (input$resultType == "num") {
      plot_cluster_sizes(expMatrix, kmeans.obj)
    }
    if (input$resultType == "expProf") {
      plot_exp_profile(expMatrix)
    }
    
    if (input$resultType == "avgExp") {
      plot_mean_exp_profile(expMatrix, kmeans.obj)
    }
    
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
