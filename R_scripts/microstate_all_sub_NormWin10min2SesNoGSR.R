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
numWinPerSub <-272
numSeed <- 4
numSub <- 21
numWinPerSeed <- numWinPerSub * numSub
totNumWin <-numWinPerSub * numSub * numSeed
# covType can be 'compCor' or 'noGSR'
covType <- 'noGSR'
winSize <- 69

sessionType <- c('2sessions')
winType <- c("Full")

		baseResults <- sprintf("/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/%s/%s", covType, sessionType)
		baseFig <- sprintf("/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/%s/%s", covType, sessionType)

winFile <- 'zWinFullCorLasso21Sub_noGSR.mat'
		fileName <- file.path(baseResults, winFile)
		matData <- readMat(fileName)
		# Copy over matrix
		varName <- 'zWinFullCorLasso21Sub'
		winAllSeeds <- matData[[varName]]
		# Check dimensions
		dim(winAllSeeds)

		# import daynamic tree cut library
		library(dynamicTreeCut)

		# clustering and cut the tree for windows of all seeds 
		distAllSeeds <- dist(winAllSeeds, method='euclidean')
		dendroAllSeeds <- hclust(distAllSeeds,method='ward')

		clustIndx <- cutreeDynamic(dendroAllSeeds, minClusterSize=2, 	distM=as.matrix(distAllSeeds))
		
		# save the index file
	indxFileName <- file.path(baseResults,sprintf("clustIndxNormWinAllSeeds_%sCorLasso_%swin%d_%s_6clusters.txt", winType, sessionType,winSize, covType))
	write.table(clustIndx,file=indxFileName,row.names=FALSE,col.names=FALSE,qmethod="double")
	

	# plot and save the dendrogram plot
		figName <- file.path(baseFig,sprintf("zWinAllSeeds_%sCorLasso_%swin%d_%s.png",winType, sessionType, winSize, covType))
		png(filename = figName , width = 480, height = 480, units = "px", pointsize = 12,bg = "white")
		figTitle <- sprintf("zWinAllSeeds_%sCorLasso_%s_%s",winType, sessionType, covType)
		plot(dendroAllSeeds,main=figTitle,ylab="Heights",xlab='FC windows', labels=F)
		dev.off()

		
		

		
		
							
		

