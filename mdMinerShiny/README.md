# Hello, this is the repository of mdMinerShiny project, we will put all code to this subfolder
# contact: robert.fh.li@gmail.com; BMI@OSU;

# How to run Shiny App
wDir <- c('/Users/li150/FHLosu/gitRepository1/')
setwd(wDir)
library(shiny)
runApp('./mdMinerShiny')

@ui-font-size: 14px **mdMiner-The Ohio State University**

mdMiner is an application for predicting drugs and drug combinations for individual cancer patients using their genomic data.
 
**Installing** 

1) Download the latest version of **RStudio** and **RShiny**.
http://www.rstudio.com/products/rstudio/download/preview/
https://www.rstudio.com/products/shiny-2/
2) Install **ShinyDashboard** package by typing the following command on the R console.  
install.packages("shinydashboard") 
3) Install networkD3 using the command 
4) Install DT by typing the following command on the R console. 
if (!require("DT")) install.packages('DT')
sessionInfo()
, source ("https://bioconductor.org/biocLite.R"); biocLite("org.Hs.eg.db”)and source("https://bioconductor.org/biocLite.R"biocLite("graphite"). 
*Instructions*
Upload fold change data of a patient in a .txt file format. Drug suggestions will be generated in descending order. Click the desired drug  in the table to display drug network and merge network. The networks display genetic and drug data. Click the title of each block to download corresponding data.
