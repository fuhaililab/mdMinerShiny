shinyUI(fluidPage(
  titlePanel("MD Miner"),
  sidebarLayout(
    sidebarPanel(
      helpText("Please choose a cancer type and a data type to download corresponding data."),

      selectInput("cancerType", "Choose A Cancer Type:", 
                  choices = c("LAML", "ACC", "BLCA", "LGG", "BRCA", "CESC", "CHOL", "COAD",
                   "ESCA", "FPPP", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", 
                   "LUSC", "DLBC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", 
                   "SKCM", "STAD", "TGCT", "THYM", "THCA", "UCS", "UCEG", "UVM")),
      selectInput("dataType", "Choose A Data Type:", 
                  choices = c("RNA-seq", "DNA methylation", "DNA copy number", 
                    "Protein expression", "miRNA-seq")),
      actionButton("action", label="Download")
    ),
    mainPanel(
      textOutput("text"),
      tableOutput("table")
    )
  )
))