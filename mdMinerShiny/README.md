# **mdMiner** 
## Repositioning Drugs for Individual Cancer Patients based on Personal Genomics using Precision Medicine

<img src="https://www.osu.edu/assets/site/images/osu-logo.png" alt="THE OHIO STATE UNIVERSITY">

### mdMiner is an innovative web application that uses machine learning approaches to predict the best drugs and drug combinations for individual cancer patients using their genomic data. mdMiner is an application developed by the Department of Biomedical Informatics at The Ohio State University College of Medicine under the supervision of Dr. Fuhai Li. 

 
### **Installing** 

1) Download the latest version of *RStudio* and *RShiny*.
https://www.rstudio.com/products/rstudio/download2/#download 
https://www.rstudio.com/products/shiny-2/

2) Install *ShinyDashboard* package by typing the following command on the R console.  
install.packages("shinydashboard") 

3) Install *RCurl* package. install.packages("RCurl")

4) Install *networkD3* https://cran.r-project.org/web/packages/networkD3/index.html

4) Install *DT* by typing the following command on the R console. 
if (!require("DT")) install.packages('DT')
sessionInfo()

5) Install *Bioconductor* package 
source ("https://bioconductor.org/biocLite.R"); biocLite("org.Hs.eg.db”)

6) Install *Graphite* package 
source("https://bioconductor.org/biocLite.R");
biocLite("graphite”)

7) Install R package *igraph* using th following command.
 install.packages("igraph")
 
### *Instructions*

Click on the *Choose File* button and upload fold change data of a patient in a .txt file format. 
Drug suggestions will be generated in descending order. 
Click the desired drug  in the table to display drug network, gene network and the merge network. 
Click the title of each block to download corresponding data.

#### For any questions, please contact: robert.fh.li@gmail.com; BMI@OSU;
