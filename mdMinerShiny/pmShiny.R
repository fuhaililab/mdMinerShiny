# for debuging: 
if (1 == 2){
	wDir <- c('/Users/li150/FHLosu/gitRepository1/mdMinerShiny/')
	setwd(wDir)
	f1 <- c("./foldchangePc3.txt")
	x1 <- read.table(f1, header=F, sep='\t')  # x1 variable has the data
	x1 <- as.matrix(x1)  # conver the list format into matrix format

	gSym <- as.character(x1[,1])  # get the gene Symbols
	fc <- as.numeric(x1[,2])  # get the value of fold change
	x = getPersonalNet1(fc, gSym);


}

# 1) install R packages
library(org.Hs.eg.db)
library(graphite)
library(igraph)

# 2) read in the fold change data
# f1 <- c("path to fold change file/file name")  # fiel name: for example: 

#fDrugMoa <- c('./drugMoaNets/')
getRepositionDrugs <- function(netPatient, n1){
	#dir1 <- c("/Users/li150/FHLosu/Projects/lincs/drugMoaNets/")
	dir1 <- c("./drugMoaNets/")
	drugList <- dir(dir1)
	drugName <- gsub('.txt', '', drugList)
	nDrug <- length(drugList)

	ds <- rep(-1.0, nDrug)
	for (i in 1:nDrug){
		netDrug <- read.table(paste(dir1, drugList[i], sep='/'), header=F, sep='\t')
		ds[i] <- getDrugRepositionScore(netPatient, netDrug)
	}
	idx <- order(ds, decreasing=T)
	ds <- ds[idx]
	drugName <- drugName[idx]

	n1 <- min(n1, nDrug)
	dat <- cbind(drugName[1:n1], ds[1:n1])

}

getDrugRepositionScore <- function(netPatient, netDrug){
	node0 <- union(netPatient[,1], netPatient[,2])
	node1 <- union(netDrug[,1], netDrug[,2])
	s <- length(intersect(node0, node1))/length(node0)

	return(s)
}


getPersonalNet2 <- function(fc, gSym, rootGenes){
	
	options(warn = -1)
	library(igraph)

	eKegg <- getKeggNet1(gSym)
	eKegg <- eKegg[,c(1,2)]  #only source/target information
	gTmp <- graph.edgelist(eKegg)  # build the background network with kegg edges

	nKegg <- union(eKegg[,1], eKegg[,2]) 

	tf1 <- getActiveTF3(fc, gSym)
	tf1 <- nKegg[nKegg %in% tf1]
	tf1 <- tf1[3]

	#get the root genes:
	# rootGenes <- read.table('./rootGenes.txt', header=F)
	# rootGenes <- as.character(rootGenes[[1]])
	rootGenes <- nKegg[nKegg %in% rootGenes]

	if (length(rootGenes) < 1){
		net1 <- matrix('test', 2,2)
		return(net1)
	}
	# get the KEGG background network
	net1 <- linkNodes1(gTmp, rootGenes, tf1)  #net1 is the network: source (1 column) and target node (2 column)

	# display the net1
	return(net1)
}


getPersonalNet1 <- function(fc, gSym){
	
	options(warn = -1)
	library(igraph)

	eKegg <- getKeggNet1(gSym)
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
	pathTmp <- pathTmp[-1,]
	return(pathTmp)
}

getKeggNet1 <- function(gSym1){
	library(org.Hs.eg.db)
	library(graphite)	
	
	entrezId <- names(as.list(org.Hs.egSYMBOL[]))
    eEntrez=lapply(kegg,function(x){return(nodes(x))})
	nSymbol=lapply(eEntrez,function(x){x=intersect(x,entrezId);unlist(as.list(org.Hs.egSYMBOL[x]))})  # node symbol

	eSymbol <- list()
	
	for (i in 1:length(kegg)){
		print(i)
		# e1 <- as.matrix(edges(kegg[[i]]))  # check the 'attributes()'
		e1 <- as.matrix(kegg[[i]]@edges)
		e1 <- e1[e1[,1] %in% entrezId & e1[,2] %in% entrezId, ]

		dim(e1) <- c(length(e1)/4, 4)

		e1[,1] <- unlist(as.list(org.Hs.egSYMBOL[e1[,1]]))
		e1[,2] <- unlist(as.list(org.Hs.egSYMBOL[e1[,2]]))
		
		e2 <- e1[e1[,3] == "undirected",]  # convert 'undirected' to directed
		if (length(e2) > 0){
			dim(e2) <- c(length(e2)/4, 4)
			e2=e2[,c(2:1,3:4)]	
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


