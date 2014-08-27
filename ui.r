


shinyUI(fluidPage(
  titlePanel("Load Estimation"),
  

  
  sidebarLayout(
    sidebarPanel( "Model Selection", 
                  radioButtons("station", "Station", 
                               choices = list("Quinapoxet", "Stillwater", "Gates")),
                  radioButtons("variable", "Variable", 
                               choices = list("NO3-N" = "NO3", "TP" = "TP", "TOC" = "TOC")),
                  dateRangeInput("dateRange", "Date Range", "2011-12-15",
                                 "2013-05-15", min = "1998-07-15", max = "2013-05-15")                  
                  ),
    
    mainPanel("Results", 
      tabsetPanel(type = "tabs", 
        tabPanel("Observations",
                 h4("Time Series"),
                 plotOutput("obsPlot"),
                 h4("Table"),
                 tableOutput("obsTable")),
        tabPanel("Load Estimates", 
           fluidRow(
             column(width = 6, 
                checkboxGroupInput("yplot1", "Show on plot", 
                   choices = list("Observations" = "meas", "Gam" = "gam", 
                      "Lm" = "lm", 
                      "Linear Interpolation - Routine samples" = "int.mon",
                      "Linear Interpolation - Routine + event samples" = "int.all",
                      "Constant Interpolation - Routine samples" = "step.mon",
                      "Constant Interpolation - Routine + event samples" = "step.all"),
                   inline = F)),
             column(width = 6, 
                radioButtons("toPlot", "What To Plot", 
                             choices = list("Concentration" = "conc", "Load" = "load")),
                checkboxInput("boundconc", "Bound predictions by observed range",
                              value = TRUE))
             
             ),
           fluidRow(
             h4("Time series"),
             plotOutput("TSplot1"),
             h3("Cumulative Load Estimates"),             
             radioButtons("seMethod", "Standard Error estimation method", 
                          choices = list("empirical" = "empirical", 
                                         "Model se.fit" = "se.fit"), inline = T),
             radioButtons("covmatQuantity", "Quantity to use for estimating covariance matrix",
                          choices = list("load" = "lhat", "concentration" = "chat"), inline = T),
             checkboxInput("addResidVar", "Include variance component from residuals", value = TRUE),
             h4("Load Table"),
             tableOutput("loadTable"),
             h4("Load Bar Plot"),
             plotOutput("loadBarPlot"))
             ),
        tabPanel("Model Info",
                 fluidRow(
                   h4("Model Forms"),
                   tableOutput("modelFormTable")
                   ),
                 fluidRow(
                   column(width = 6,
                          radioButtons("termVarLm", "LM variable to inspect",
                          choices = c("label 1" = "option1", "label 2" = "option2"))
                          ),
                   column(width = 6, 
                          radioButtons("termVarGam", "GAM variable to inspect",
                            choices = c("label 1"= "option1", "label 2" = "option2"))
                          )
                   ),
                 fluidRow(
                   column(width = 6,
                          h4("LM Term Plots"), 
                          plotOutput("termPlotsLM")),
                   column(width = 6, 
                          h4("GAM Term Plots"),
                          plotOutput("termPlotsGAM"))
                   ),
                 fluidRow(
                   h3("Fit Diagnostics"),
                   plotOutput("fitPlotConc"),
                   plotOutput("fitPlotLoad")
                   ),
                 fluidRow(h3("Residual Diagnostics")),
                 fluidRow(
                   column(width = 6, 
                    radioButtons("residXvarLm", "LM variable for x-axis",
                      choices = c("label 1" = "option1", "label 2" = "option2"))
                   ),
                   column(width = 6, 
                    radioButtons("residXvarGam", "GAM variable for x-axis",
                      choices = c("label 1" = "option1", "label 2" = "option2"))
                   )),
                 fluidRow(
                   column(width = 6,
                    h4("LM Residual Plots"),
                    plotOutput("lmResidPlot"),
                    plotOutput("qqPlotLM"),
                    plotOutput("boxplotLM"),
                    plotOutput("acfLM")),
                   column(width = 6,
                    h4("GAM Residual Plots"),
                    plotOutput("gamResidPlot"),
                    plotOutput("qqPlotGAM"),
                    plotOutput("boxplotGAM"),
                    plotOutput("acfGAM"))
                 )
                 ),
        tabPanel("Results Download",
          downloadButton('downloadData', 'Download'),
          tableOutput('resultsTS')),
        tabPanel("Test",
#           tableOutput("test"),
          textOutput("test2"))
#           tableOutput("test2"),
#           plotOutput("test3"))
        )
      )
    )
  )
)