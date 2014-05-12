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
sessionType <- 'session1'
covType <- 'GSR'
baseResults <- sprintf("/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/%s/%s/21seedTS/", covType, sessionType)
baseFig <- sprintf("/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/%s/%s/21seedTS/", covType, sessionType)
		fileName <- file.path(baseResults, 'normMeanTS21seed.mat')
		matData <- readMat(fileName)
		# Copy over matrix
		varName <- 'normMean2D'
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
	indxFileName <- file.path(baseResults,sprintf("clustIndx21seedTS_%s_%s.txt", sessionType,covType))
	write.table(clustIndx,file=indxFileName,row.names=FALSE,col.names=FALSE,qmethod="double")
	
	# plot and save the dendrogram plot
		figName <- file.path(baseFig,sprintf("zMeanTS21Seed_%s_%s.png",sessionType, covType))
		png(filename = figName , width = 480, height = 480, units = "px", pointsize = 12,bg = "white")
		figTitle <- sprintf("zMeanTS21Seed_%s_%s",sessionType, covType)
		plot(dendroAllSeeds,main=figTitle,ylab="Heights",xlab='21 seeds')
		dev.off()
		

		
		
							
		

