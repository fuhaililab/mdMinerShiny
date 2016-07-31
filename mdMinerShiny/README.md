## <img src="https://www.osu.edu/assets/site/images/osu-logo.png" alt="THE OHIO STATE UNIVERSITY">
## Repositioning Drugs for Individual Cancer Patients based on Personal Genomics using Precision Medicine

### mdMiner is an innovative web application that uses machine learning approaches to predict the best drugs and drug combinations for individual cancer patients using their genomic data. mdMiner is an application developed by the Department of Biomedical Informatics at The Ohio State University College of Medicine under the supervision of Dr. Fuhai Li. 

 
### **Installing**  
1) Install *R* on your machine following the instructions on any of these links.  
  Windows - https://cran.r-project.org/bin/windows/base/  
  Linux - https://cran.r-project.org/doc/manuals/r-release/R-admin.html  
  Mac OS X - https://cran.r-project.org/bin/macosx/

2) Download and install the latest version of *RStudio*  
   https://www.rstudio.com/products/rstudio/download2/#download 

3) Install *Shiny* and *ShinyDashboard* packages by typing the following commands on the R console.  
  \> install.packages("shiny");  
  \> install.packages("shinydashboard");  
 
4) Install other dependencies by typing the following commands on the R console.  
  \> install.packages("RCurl");  
  \> install.packages("networkD3");  
  \> install.packages("httr");  
  \> install.packages("DT");  
  \> install.packages("igraph");

5) Use *Bioconductor* as following to install *Graphite* and *Genome annotation* packages. Type the following commands in R console.   
  \> source("https://bioconductor.org/biocLite.R");  
  \> biocLite();  
  \> biocLite("graphite”);  
  \> biocLite("org.Hs.eg.db”);  

### **Instructions**  
Click *Run App* in RStudio. This will open the webapp in a browser.  
Click on the *Choose File* button and upload fold change data of a patient in a .txt file format.  
Drug suggestions will be generated in descending order.  
Click the desired drug  in the table to display drug network, gene network and the merge network.  
Click the title of each block to download corresponding data.  

#### For any questions, please contact: robert.fh.li@gmail.com.
#####- Dr. Fuhai Li - Department of Biomedical Informatics, The Ohio State University.

