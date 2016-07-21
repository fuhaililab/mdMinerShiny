# Hello, this is the repository of mdMinerShiny project, we will put all code to this subfolder
# contact: robert.fh.li@gmail.com; BMI@OSU;

# How to run Shiny App
wDir <- c('/Users/li150/FHLosu/gitRepository1/')
setwd(wDir)
library(shiny)
runApp('./mdMinerShiny')

mdMiner is an application for predicting drugs and drug combinations for individual cancer patients using their genomic data. 
**Installing**
Install Shiny, ShinyDashboard, DT, source ("https://bioconductor.org/biocLite.R"); biocLite("org.Hs.eg.db”)and source("https://bioconductor.org/biocLite.R"biocLite("graphite"). 
*Instructions*
Upload fold change data of a patient in a .txt file format. Drug suggestions will be generated in descending order. Click the desired drug  in the table to display drug network and merge network. The networks display genetic and drug data. Click the title of each block to download corresponding data. 
