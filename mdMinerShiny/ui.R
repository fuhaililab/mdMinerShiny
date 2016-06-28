library(networkD3);
data(MisLinks)
data(MisNodes)
shinyUI(navbarPage("MD-Miner Shiny",

  tabPanel("Download Data",
    titlePanel("Download Data You Need from TCGA"),
    sidebarPanel(
      helpText("Please choose the type of cancer and data to download."),
      selectizeInput("cancerType", "Choose Cancer Types:", 
                  choices = c("LAML", "ACC", "BLCA", "LGG", "BRCA", "CESC", "CHOL", "COAD",
                   "ESCA", "FPPP", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", 
                   "LUSC", "DLBC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", 
                   "SKCM", "STAD", "TGCT", "THYM", "THCA", "UCS", "UCEG", "UVM"), multiple = TRUE),
      selectInput("dataType", "Choose A Data Type:", 
                  choices = c("RNA-seq", "DNA methylation", "DNA copy number", 
                    "Protein expression", "miRNA-seq", "Clinical"), multiple = TRUE),
      downloadButton("downloadData", label="Download"),
      helpText("Data will be stored in your /Downloads file")
    )
  ),

  tabPanel("Comprehensive Analysis",
    titlePanel("Analyze patient data with existing samples"),

    sidebarLayout(
      sidebarPanel(
        fileInput('file1', 'Choose Patient Fold Change Data File', 
          accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
      ),
    
      mainPanel(
        tabsetPanel(
          tabPanel("Network Display", forceNetworkOutput("force")), 
          tabPanel("Cluster Analysis", plotOutput("plot")), 
          tabPanel("Survival Prediction", tableOutput("table"))
        )
      )
    )
  ),

  tabPanel("Drug Suggestion", 
    titlePanel("Predict the best drug according to genetic data")
  )
))
