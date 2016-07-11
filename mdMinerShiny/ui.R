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
                  accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
        tags$h5("Click to Download the Top Drug Suggestions Table"),
        downloadButton('downloadData', 'Download')
      ),
      dashboardBody( 
        fluidRow(
          box(title=tags$b("Top Drug Suggestions"), value = tags$p(style = "font-size: 10px;", tags$b()),  solidHeader= TRUE,collapsible = TRUE,  status = "info", DT::dataTableOutput("table", width = 924, height = 200), width =12)
          ),
        fixedRow(
          box(title=tags$b("Patient and Drug Merged Netowrk"), value = tags$p(style = "font-size: 10px;"), solidHeader= TRUE, status = "primary", collapsible = TRUE,width =12, verbatimTextOutput("text"), forceNetworkOutput("mergeNetwork", height = 400, width = 900))
        ),
        fluidRow(
          box(title=tags$b("Patient Information and Gene Network"), value = tags$p(style = "font-size: 10px;"), solidHeader= TRUE, status = "primary", collapsible = TRUE, forceNetworkOutput("patientNetwork", height = 400, width = 400)),
          box(title=tags$b("Drug Suggestion and Gene Network"), value = tags$p(style = "font-size: 10px;"), solidHeader= TRUE, status = "primary", collapsible = TRUE, forceNetworkOutput("drugNetwork", height = 400, width = 400))
          )
        # fixedRow(
        #   column(6,
        #     forceNetworkOutput("force")
        #   ),
        #   column(6,
        #     forceNetworkOutput("geneforce")
        #   ),
          )
        )
      )
    )
  
