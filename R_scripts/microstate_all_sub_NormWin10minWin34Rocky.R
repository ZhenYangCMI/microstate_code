#clear screen
command+alt+l

#clear workspace
rm(list = ls(all = TRUE))

###
# Load Data
###

# Load matlab library
library(R.matlab)

# Read in matlab data
# numWinPerSub 272 for winSize69, 284 for winSize 34, and 250 for winSize 136
# they corresponding to winSize 44s, 22s, and 88s
winSize <- 136
numWinPerSub <-284
numSeed <- 4
numSub <- 22
numWinPerSeed <- numWinPerSub * numSub
totNumWin <-numWinPerSub * numSub * numSeed

sessionType <- c('session1')
sessionNum <- 1:length(sessionType)
winType <- c("Full")
winTypeNum <- 1:length(winType)
seed <- c('seed1','seed2','seed3','seed4')
seedNum <- 1:length(seed)


for (i in sessionNum) {
	for (j in winTypeNum) {
		baseResults <- sprintf("/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/%s", sessionType[i])
		baseFig <- sprintf("/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/%s", sessionType[i])

		#fileName <- file.path(baseResults, "test.mat")
		fileName <- file.path(baseResults, sprintf("zWin%sCorLasso_OptimalLambdaPerSub_645_%swin%d.mat", winType[j], sessionType[i], winSize))
		matData <- readMat(fileName)
		# Copy over matrix
		varName <- sprintf("zWin%sCorLasso",winType[j])
		winAllSeeds <- matData[[varName]]
		# Check dimensions
		dim(winAllSeeds)

		# import daynamic tree cut library
		library(dynamicTreeCut)

		# clustering and cut the tree for windows of all seeds 
		distAllSeeds <- dist(winAllSeeds, method='euclidean')
		dendroAllSeeds <- hclust(distAllSeeds,method='ward')
		clustIndx <- cutreeDynamic(dendroAllSeeds, minClusterSize=2,	distM=as.matrix(distAllSeeds))
		
		# save the index file
	indxFileName <- file.path(baseResults,sprintf("clustIndxWinAllSeeds_%sCorLasso_%s_10min_win%d.txt", winType[j], sessionType[i], winSize)	)
	write.table(clustIndx,file=indxFileName,row.names=FALSE,col.names=FALSE,qmethod="double")
	
	# plot and save the dendrogram plot
		figName <- file.path(baseFig,sprintf("zWinAllSeeds_%sCorLasso_%s_10min_win%d.png",winType[j], sessionType[i], winSize))
		png(filename = figName , width = 480, height = 480, units = "px", pointsize = 12,bg = "white")
		figTitle <- sprintf("zWinAllSeeds_%sCorLasso_%s_10min_win%d",winType[j], sessionType[i], winSize)
		plot(dendroAllSeeds,main=figTitle,ylab="Heights",xlab='FC windows', labels=F)
		dev.off()
		
	}
	}




