library(shiny)
library(networkD3)
library(shinydashboard)
library(DT)

shinyUI(fluidPage(theme = "bootstrap.css", 
    tags$style(type="text/css",
      "label {font-size: 12px;}",
      ".recalculating {opacity: 1.0;}"),
   
    
    
    dashboardPage(skin = "red",
      dashboardHeader(
        title = tags$h6("MdMiner: The Ohio State University"),
        titleWidth= 230
      ),
      dashboardSidebar(
        #width = 450,
        fileInput('file1',tags$h5("Choose Patient Fold Change Data"), 
                  accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')) #,
        #tags$h5("Top Drug Suggestions"),
        #DT::dataTableOutput("table")
      ),
      dashboardBody( 
        fluidRow(
          box(title= tags$h5("Top Drug Suggestions"), solidHeader= TRUE, status = "info", collapsible = TRUE, DT::dataTableOutput("table", width = 924, height = 200), width =12)
          ),
        fluidRow(
          box(title= tags$h5("Patient Information and Gene Network"), solidHeader= TRUE, status = "primary", collapsible = TRUE, forceNetworkOutput("force", height = 400, width = 400)),
          box(title =tags$h5("Drug Suggestion and Gene Network"), solidHeader= TRUE, status = "primary", collapsible = TRUE) #, verbatimTextOutput("text"), height = 400)
          ),
        # fixedRow(
        #   column(6,
        #     forceNetworkOutput("force")
        #   ),
        #   column(6,
        #     forceNetworkOutput("geneforce")
        #   ),
        fixedRow(
          box(title= tags$h5("Patient and Drug Merged Netowrk"), solidHeader= TRUE, status = "primary", collapsible = TRUE),
          box(title= tags$h5("Survival Analysis"), solidHeader= TRUE, status = "primary", collapsible = TRUE)
        ),
        fixedRow(
          column(6,
                 forceNetworkOutput("mergeforce")
          ),
          column(6,
                 plotOutput("survival")
          )
        )
      )
    )
  )
)