# 1) install R packages
# library(org.Hs.eg.db)
# library(graphite)
# library(igraph)

# 2) read in the fold change data
# f1 <- c("path to fold change file/file name")  # fiel name: for example: f1 <- c("./foldchangePc3.txt")
# x1 <- read.table(f1, header=F, sep='\t')  # x1 variable has the data
# x1 <- as.matrix(x1)  # conver the list format into matrix format

# gSym <- as.character(x1[,1])  # get the gene Symbols
# fc <- as.numeric(x1[,2])  # get the value of fold change

getPersonalNet1 <- function(fc, gSym){
	library(igraph)

	eKegg <- getKeggNet1()
	eKegg <- eKegg[,c(1,2)]  #only source/target information
	gTmp <- graph.edgelist(eKegg)  # build the background network with kegg edges

	nKegg <- union(eKegg[,1], eKegg[,2]) 

	tf1 <- getActiveTF3(fc, gSym)
	tf1 <- nKegg[nKegg %in% tf1]
	tf1 <- tf1[3]


	#get the root genes:
	rootGenes <- read.table('./rootGenes.txt', header=F)
	rootGenes <- as.character(rootGenes[[1]])
	rootGenes <- nKegg[nKegg %in% rootGenes]

	# get the KEGG background network
	net1 <- linkNodes1(gTmp, rootGenes, tf1)  #net1 is the network: source (1 column) and target node (2 column)

	# display the net1
	return(net1)
}

linkNodes1 <- function(gTmp, recTmp, tfTmp){
	library(igraph)
	
	vTmp <- V(gTmp)$name
	nPath1 <- 0
	pathTmp <- matrix('test', 1,2)
	vt <- matrix('test', 1,2)

	nRecTmp <- length(recTmp)
	for (j in 1:nRecTmp){
		paths <- get.shortest.paths(gTmp, recTmp[j], tfTmp, mode='out')
		paths <- paths$vpath 
		nPath <- length(paths)
		
		
		if (nPath > 0){
			for (k in 1:nPath){
				pt <- paths[[k]]
				pt <- pt$name
				nPt <- length(pt)
				if (nPt < 2){
					next	
				}
				
				for (l in 1:(nPt-1)){
					#vt <- vTmp[c(pt[l], pt[l+1])]
					vt <- c(pt[l], pt[l+1])
					pathTmp <- rbind(pathTmp, vt) 	
				}
			}
		}
	}
	return(pathTmp)
}

getKeggNet1 <- function(){
	library(org.Hs.eg.db)
	library(graphite)	
	
	entrezId <- names(as.list(org.Hs.egSYMBOL[]))
	eEntrez=lapply(kegg,function(x){return(nodes(x))})	
	nSymbol=lapply(eEntrez,function(x){x=intersect(x,gSym1);unlist(as.list(org.Hs.egSYMBOL[x]))})  # node symbol

	eSymbol <- list()
	
	for (i in 1:length(kegg)){
		e1 <- as.matrix(edges(kegg[[i]]))
		e1 <- e1[e1[,1] %in% entrezId & e1[,2] %in% entrezId, ]

		dim(e1) <- c(length(e1)/4, 4)

		e1[,1] <- unlist(as.list(org.Hs.egSYMBOL[e1[,1]]))
		e1[,2] <- unlist(as.list(org.Hs.egSYMBOL[e1[,2]]))
		
		e2 <- e1[e1[,3] == "undirected",]  # convert 'undirected' to directed
		if (length(e2) > 0){
			dim(y) <- c(length(y)/4, 4)
			y=y[,c(2:1,3:4)]	
			e1 <- rbind(e1, e2)
		}
		eSymbol[[i]] <- e1[!duplicated(e1),]
	}
	
	names(eSymbol)=names(kegg)

	n1 <- length(eSymbol)
	e1 <- eSymbol[[1]]
	for (i in 2:n1){
		e1 <- rbind(e1, eSymbol[[i]])
	}

	return(e1)

}

getActiveTF3 <- function(fc, gSym){
	
	tfTar0 <- read.table('./tfTarget.txt', header=F, sep='\t')
	tfTar0 <- as.matrix(tfTar0)

	# identify the activated TFs
	tar1 <- unique(tfTar0[,2])
	idxt1 <- which(gSym %in% tar1)
	
	sTar1 <- fc[idxt1] # suppression (negative zscore) for drugs
	
	tf1 <- unique(tfTar0[,1])
	nTf <- length(tf1)
	
	sTf <- rep(0, nTf)  # importance score of TF
	for (i in 1:nTf){
		str1 <- tf1[i]  # i-th TF
		str2 <- tfTar0[tfTar0[,1] == str1, 2]  #targets of given TF
		idxt2 <- which(tar1 %in% str2)
		st <- sTar1[idxt2]
		st1 <- st[order(st, decreasing=T)]
		nt1 <- min(3, length(st1))
		st1 <- st1[1:nt1]
		sTf[i] <- sum(st1)/nt1	
	}
	
	idxt1 <- order(sTf, decreasing=T)
	sTf <- sTf[idxt1]
	tf2 <- tf1[idxt1]


	tf2 <- tf2[1:15]
	
	return(tf2)
}


# f1 <- paste(dirProject1, '/PC3/foldchangePc3.txt', sep='')
# write.table(fcPc3, f1, row.names=T, col.names=F, sep='\t')
# f1 <- paste(dirProject1, '/PC3/tfTarget.txt', sep='')
# write.table(tfTar0, f1, row.names=F, col.names=F, sep='\t')
# tfTar0 <- getTfTargetInteraction1(dirDat0)


