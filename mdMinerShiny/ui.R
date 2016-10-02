library(shiny)
library(networkD3)
library(shinydashboard)
library(DT)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("About", tabName = "about", icon = icon("dashboard")),

    menuItem("Start", icon = icon("th"), tabName = "start",
      fileInput('file1',tags$h5(tags$b("Choose Patient Fold Change Data")), 
                accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
    ), 
    
    menuItem("Drug ranking", tabName = "table", icon = icon("dashboard")),

    menuItem("Network Display", icon = icon("th"), tabName = "network",
    	menuItem("Patient Network", icon = icon("th"), tabName = "network1"),
    	menuItem("Drug Network", icon = icon("th"), tabName = "network2"),
    	menuItem("Combined Network", icon = icon("th"), tabName = "network3")
    )
))

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "about",
      h2("Detailed description about MD-Miner"),
      h2("Instructions")
    ),

    tabItem(tabName = "table",
      fluidRow(
        box(title=downloadButton('downloadData', tags$b("Top Drug Suggestions")), value = tags$p(style = "font-size: 10px;", tags$b()),  solidHeader= TRUE,collapsible = TRUE, DT::dataTableOutput("table"), width =12)
        )
    ),

    tabItem(tabName = "network1",
      fluidRow(
        box(width = 12, title=downloadButton('downloadData1', tags$b("Patient Information and Gene Network")), value = tags$p(style = "font-size: 10px;"), solidHeader= TRUE, collapsible = TRUE, forceNetworkOutput("patientNetwork", width = 800))
        )
    ),

    tabItem(tabName = "network2",
      fluidRow(
        box(width = 12, title=downloadButton('downloadData2', tags$b("Drug Network")), value = tags$p(style = "font-size: 10px;"), solidHeader= TRUE, collapsible = TRUE,  forceNetworkOutput("drugNetwork", width = 800))
        )
    ),

    tabItem(tabName = "network3",
      fluidRow(
        box(width = 12, title=downloadButton('downloadData3', tags$b("Drug Suggestion and Gene Network")), value = tags$p(style = "font-size: 10px;"), solidHeader= TRUE, collapsible = TRUE, forceNetworkOutput("mergeNetwork", width = 800))
        )
    )

  )
)

shinyUI(fluidPage(theme = "bootstrap.css", 
    tags$style(type="text/css",
      "label {font-size: 12px;}",
      #"label {font-weight: bold;}",
      ".recalculating {opacity: 1.0;}"
    ),
    
    dashboardPage(skin = "yellow",
      dashboardHeader(title = tags$b("MD-Miner")),
	    sidebar,
	    body
    )

))