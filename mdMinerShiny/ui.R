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
        tags$h5("App Description")
      ),
      dashboardBody( 
        fluidRow(
          box(title=downloadButton('downloadData', tags$b("Top Drug Suggestions")), value = tags$p(style = "font-size: 10px;", tags$b()),  solidHeader= TRUE,collapsible = TRUE,  status = "info", DT::dataTableOutput("table"), width =12)
          ),
        fixedRow(
          box(title=downloadButton('downloadData1', tags$b("Patient and Drug Merged Netowrk")), value = tags$p(style = "font-size: 10px;"), solidHeader= TRUE, status = "primary", collapsible = TRUE,width =12, verbatimTextOutput("text"), forceNetworkOutput("mergeNetwork"))
        ),
        fluidRow(
          box(title=downloadButton('downloadData2', tags$b("Patient Information and Gene Network")), value = tags$p(style = "font-size: 10px;"), solidHeader= TRUE, collapsible = TRUE, status = "primary", forceNetworkOutput("patientNetwork")),
          box(title=downloadButton('downloadData3', tags$b("Drug Suggestion and Gene Network")), value = tags$p(style = "font-size: 10px;"), solidHeader= TRUE, status = "primary", collapsible = TRUE, forceNetworkOutput("drugNetwork"))
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
  
