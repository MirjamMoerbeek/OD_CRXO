library(shinythemes)
library(plotly)
library(shinydashboard)
library(DT)

navbarPage(theme = shinytheme("cerulean"),"Optimal design of cluster randomized crossover trials with time periods of fixed duration ",
           tabPanel("Optimal number of time periods",
                    sidebarLayout(
                      sidebarPanel(
                                   h4("Correlation parameters"),
                                   sliderInput("ICC1", "Base correlation", min = 0, max = 0.2, value = 0.036, step = 0.001),
                                   sliderInput("CAC1", "Cluster autocorrelation", min = 0, max = 1, value = 0.77^(1/11), step = 0.01),
                                   
                                   h4("Trial configuration"),
                                   sliderInput("nr.subj1", "Number of subjects per cluster-period", min = 0, max = 500, value = 50, step = 1),
                                   sliderInput("nr.per1", "Minimum and maximum number of time periods:", min = 2, max = 12, value = c(2,12),step=1),
                                   
                                   h4("Costs and budget"),
                                   numericInput("costs.c1", label = h5("Costs per cluster"), value = 2500),
                                   numericInput("costs.s1", label = h5("Costs per subject"), value = 50),
                                   numericInput("costs.x1", label = h5("Costs per teatment switch"), value = 250),
                                   numericInput("budget1", label = h5("Budget"), value = 2500000),
                                   
                                   submitButton("Submit")
                                  ),
                      
                      mainPanel(h3(""),
                                tabBox(
                                  title = "",width=8,  height = 550,
                                  # The id lets us use input$tabset1 on the server to find the current tab
                                  id = "tabset1", 
                                  tabPanel("Graph", 
                                           plotlyOutput("plot1", width = 800, height = 700),
                                  ),
                                  
                                  tabPanel("Table", 
                                          DT::dataTableOutput(outputId = "ResultsTable")       
                                  )
                                )
                      )
                    )
           ),
           
           tabPanel("Optimal number of treatment switches",
                    sidebarLayout(
                      sidebarPanel(
                        h4("Correlation parameters"),
                        sliderInput("ICC2", "Base correlation", min = 0, max = 0.2, value = 0.036, step = 0.001),
                        sliderInput("CAC2", "Cluster autocorrelation", min = 0, max = 1, value = 0.77^(1/11), step = 0.01),
                        
                        h4("Trial configuration"),
                        sliderInput("nr.subj2", "Number of subjects per cluster-period", min = 0, max = 500, value = 50, step = 1),
                        sliderInput("nr.per2", "Number of time periods:", min = 2, max = 12, value = 4, step=1),
                        
                        h4("Costs and budget"),
                        numericInput("costs.c2", label = h5("Costs per cluster"), value = 2500),
                        numericInput("costs.s2", label = h5("Costs per subject"), value = 50),
                        numericInput("costs.x2", label = h5("Costs per teatment switch"), value = 250),
                        numericInput("budget2", label = h5("Budget"), value = 2500000),
                        
                        submitButton("Submit")
                      ),
                      
                      mainPanel(h3(""),
                                   DT::dataTableOutput(outputId = "ResultsTable2")   
                      )
                    )
           ),
     )
