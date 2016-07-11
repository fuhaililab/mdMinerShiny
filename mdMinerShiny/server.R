library(shiny)
library(networkD3)
library(DT)

shiny.maxRequestSize=30*1024^2

source("./Module_A.R");
source("./pmshiny.R");
shinyServer(function(input, output) {
	# RNASeq <- reactive({
	# 	result = NULL;
	# 	parentURL = "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/";
	# 	URL = paste(parentURL, tolower(input$cancerType), sep = "");
	# 	if (input$dataType == "RNA-seq") {
	# 		URL = paste(URL, "/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/", sep = "");
	# 		URLcontent = getTCGA_URL(URL);
	# 		dirURL = URLcontent$dir_url;
	# 		level_3_URL = grep("Level_3", dirURL, value = TRUE);
	# 		subURLcontent = getTCGA_URL(level_3_URL[length(level_3_URL)]);
	# 		fileURL = subURLcontent$file_url;
	# 		normalizedResultsURL = grep("rsem.genes.normalized_results", fileURL, value = TRUE);

	# 		withProgress(message = 'Collecting samples', {
	# 			n = 5;
	# 			# n = length(normalizedResultsURL);
	# 			for (j in 1:n) {
	# 				incProgress(1/n, detail = paste(j, '/', n))
	# 				x = read.table(normalizedResultsURL[j], skip=1);
	# 				names(x)= c("Gene ID", "patient 1");
	# 				if (is.null(result)) {
	# 					result = x;
	# 				} else {
	# 					result = cbind(result, x[, 2]);
	# 				}
	# 			}
	# 		})

	# 		return(cor(result[sapply(result, is.numeric)]));			
	# 	}

	# })

	# output$text <- renderText({
	# 	paste("You are choosing the", input$dataType, "of", input$cancerType)
	# })

	# output$downloadData <- downloadHandler(
	# 	filename = function() {
	# 		paste(input$cancerType, '_', input$dataType, '.txt', sep='')
	# 	},
	# 	content = function(file) {
	# 		write.table(RNASeq(), file, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE)
	# 	}
	# )

	networkAndDrugScore <- reactive({
	    inFile <- input$file1
	    if (!is.null(inFile)) {	    	
			foldChangePC <- read.table(inFile$datapath)
			foldChangePC <- as.matrix(foldChangePC)
			gSym <- as.character(foldChangePC[, 1])	
			fc <- as.numeric(foldChangePC[, 2])
			x <- getPersonalNet1(fc, gSym)
			y <- getRepositionDrugs(x, 20)
			Drug <- y[, 1]
			Score <- y[, 2]
			drugAndScore <- data.frame(Drug, Score)
			return(list(network = x, drugAndScore = drugAndScore))
		}
	})
	
	output$patientNetwork <- renderForceNetwork({
	  #   inFile <- input$file1;
	  #   if (!is.null(inFile)) {	    	
			# foldChangePC = read.table(inFile$datapath);
			# foldChangePC = as.matrix(foldChangePC)
			# gSym = as.character(foldChangePC[, 1])
			# fc = as.numeric(foldChangePC[, 2])
			# x = getPersonalNet1(fc, gSym);
			if (!is.null(networkAndDrugScore())) {
				x = networkAndDrugScore()$network;
				sourceName = x[, 1];
				targetName = x[, 2];
				name = unique(c(as.character(sourceName), as.character(targetName)));
				group = numeric(length(name)) + 1;
				# group[match(('CREBBP'), name)] = 2;
				# group[match(('ATF4'), name)] = 3;

				size = numeric(length(name));
				# size[match(('ATF4'), name)] = 15;
				MisNodes = data.frame(name, group, size);

				source = c(match(sourceName[1], name) - 1);
				target = c(match(targetName[1], name) - 1);
				for (i in 2:length(sourceName)) {
					source = c(source, match(sourceName[i], name) - 1);
					target = c(target, match(targetName[i], name) - 1);
				}

				value = numeric(length(source)) + 1;
				MisLinks = data.frame(source, target, value);

				forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source", Target = "target",
				            Value = "value", NodeID = "name", Nodesize = "size", Group = "group", colourScale = JS("d3.scale.category10()"),
				            linkDistance = 100, opacity = 1, opacityNoHover = 1, charge = -100, zoom = TRUE);			
			}    	
	})

	output$drugNetwork <- renderForceNetwork({
		if (!is.null(input$table_rows_selected)) {
			drugName = networkAndDrugScore()$drugAndScore[as.numeric(input$table_rows_selected), 1];
			print(drugName)
			x = read.table(paste("./drugMoaNets/", as.character(drugName), '.txt', sep=''));

			sourceName = x[, 1];
			targetName = x[, 2];
			name = unique(c(as.character(sourceName), as.character(targetName)));
			group = numeric(length(name)) + 1;
			# group[match(('CREBBP'), name)] = 2;
			# group[match(('ATF4'), name)] = 3;

			size = numeric(length(name));
			# size[match(('ATF4'), name)] = 15;
			MisNodes = data.frame(name, group, size);

			source = c(match(sourceName[1], name) - 1);
			target = c(match(targetName[1], name) - 1);
			for (i in 2:length(sourceName)) {
				source = c(source, match(sourceName[i], name) - 1);
				target = c(target, match(targetName[i], name) - 1);
			}

			value = numeric(length(source)) + 1;
			MisLinks = data.frame(source, target, value);

			forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source", Target = "target",
			            Value = "value", NodeID = "name", Nodesize = "size", Group = "group", colourScale = JS("d3.scale.category10()"),
			            linkDistance = 100, opacity = 1, opacityNoHover = 1, charge = -100, zoom = TRUE);			
		}
	})
	
	output$mergeNetwork <- renderForceNetwork({
		if (!is.null(networkAndDrugScore())) {
			x = networkAndDrugScore()$network;
			sourceName = x[, 1];
			targetName = x[, 2];
			name = unique(c(as.character(sourceName), as.character(targetName)));
			size = numeric(length(name));
			# size[match(('ATF4'), name)] = 15;
			group = numeric(length(name)) + 1;
			# group[match(('CREBBP'), name)] = 2;
			# group[match(('ATF4'), name)] = 3;
			if (!is.null(input$table_rows_selected)) {
				drugName = networkAndDrugScore()$drugAndScore[as.numeric(input$table_rows_selected), 1];
				drugNetwork = read.table(paste("./drugMoaNets/", as.character(drugName), '.txt', sep=''));
				sourceNameInDrugNetwork = drugNetwork[, 1];
				targetNameInDrugNetwork = drugNetwork[, 2];
				nameInDrugNetwork = unique(c(as.character(sourceNameInDrugNetwork), as.character(targetNameInDrugNetwork)));
				overlap = intersect(name, nameInDrugNetwork);
				print(length(overlap));
				if (length(overlap) > 0) {
					for (i in 1:length(overlap)) {
						group[match((as.character(overlap[i])), name)] = 2;
					}						
				}
				MisNodes = data.frame(name, group, size);

				source = c(match(sourceName[1], name) - 1);
				target = c(match(targetName[1], name) - 1);
				for (i in 2:length(sourceName)) {
					source = c(source, match(sourceName[i], name) - 1);
					target = c(target, match(targetName[i], name) - 1);
				}

				value = numeric(length(source)) + 1;
				MisLinks = data.frame(source, target, value);

				forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source", Target = "target",
				            Value = "value", NodeID = "name", Nodesize = "size", Group = "group", colourScale = JS("d3.scale.category10()"),
				            linkDistance = 100, opacity = 1, opacityNoHover = 1, charge = -100, zoom = TRUE);				
			}		
		}	    	

	})
	
	output$table <- DT::renderDataTable(networkAndDrugScore()$drugAndScore, server = FALSE, selection = 'single', 
	         options = list(pageLength = 5))
	
	output$downloadData <-  downloadHandler(
		filename = function() {paste("Top Drugs", ".csv", sep="")},
		content = function(file) {
			write.csv(networkAndDrugScore()$drugAndScore, file, row.names=F)
		}
	)
	
	# output$downloadData <- downloadHandler(
	# 	filename = function() {
	# 		paste(input$cancerType, '_', input$dataType, '.txt', sep='')
	# 	},
	# 	content = function(file) {
	# 		write.table(RNASeq(), file, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE)
	# 	}
	# )

	output$text <- renderPrint(
		paste('You Choose Drug: ', networkAndDrugScore()$drugAndScore[as.numeric(input$table_rows_selected), 1], sep='')
	)
})




