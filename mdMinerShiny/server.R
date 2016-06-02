library(shiny);
source("./Module_A.r");
shinyServer(function(input, output) {
	RNASeq <- reactive({
		result = NULL;
		parentURL = "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/";
		URL = paste(parentURL, tolower(input$cancerType), sep = "");
		if (input$dataType == "RNA-seq") {
			URL = paste(URL, "/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/", sep = "");
			URLcontent = getTCGA_URL(URL);
			dirURL = URLcontent$dir_url;
			level_3_URL = grep("Level_3", dirURL, value = TRUE);
			subURLcontent = getTCGA_URL(level_3_URL[length(level_3_URL)]);
			fileURL = subURLcontent$file_url;
			normalizedResultsURL = grep("rsem.genes.normalized_results", fileURL, value = TRUE);

			withProgress(message = 'Collecting samples', {
				n = 50;
				# n = length(normalizedResultsURL);
				for (j in 1:n) {
					incProgress(1/n, detail = paste(j, '/', n))
					x = read.table(normalizedResultsURL[j]);
					if (is.null(result)) {
						result = x;
					} else {
						result = cbind(result, x[, 2]);
					}
				}
			})

			return(result);				
		}

	})

	output$text <- renderText({
		paste("You are choosing the", input$dataType, "of", input$cancerType)
	})

	output$downloadData <- downloadHandler(
		filename = function() {
			paste(input$cancerType, '_', input$dataType, '.txt', sep='')
		},
		content = function(file) {
			write.table(RNASeq(), file, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE)
		}
	)

})