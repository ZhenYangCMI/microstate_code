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
base <- "/Users/zhen.yang/Documents/Zhen_CMI/microstate/results"
numWinPerSub <-272
numSeed <- 4
numSub <- 22
numWinPerSeed <- numWinPerSub * numSub
totNumWin <-numWinPerSub * numSub * numSeed

sessionType <- c('session1','session2')
sessionNum <- 1:length(sessionType)
winType <- c("Full")
winTypeNum <- 1:length(winType)
seed <- c('seed1','seed2','seed3','seed4')
seedNum <- 1:length(seed)

for (i in sessionNum) {
	for (j in winTypeNum) {
		fileName <- file.path(base, sprintf("zWin%sCorLasso_OptimalLambdaPerSub_645_session%i_10min.mat", winType[j], i))
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
		clustIndx <- cutreeDynamic(dendroAllSeeds, minClusterSize=2, 	distM=as.matrix(distAllSeeds))
		
		# save the index file
	indxFileName <- sprintf("/Users/zhen.yang/Documents/Zhen_CMI/microstate/results/lambdaOptimalPerSub/zWinAllSeeds_%sCorLasso_%s_10min.txt", winType[j], sessionType[i])	
	write.table(clustIndx,file=indxFileName,row.names=FALSE,col.names=FALSE,qmethod="double")
	
	# plot and save the dendrogram plot
		figName <- sprintf("/Users/zhen.yang/Documents/Zhen_CMI/microstate/fig/lambdaOptimalPerSub/zWinAllSeeds_%sCorLasso_%s_10min.png",winType[j], sessionType[i])
		png(filename = figName , width = 480, height = 480, units = "px", pointsize = 12,bg = "white")
		figTitle <- sprintf("zWinAllSeeds_%sCorLasso_%s_10min",winType[j], sessionType[i])
		plot(dendroAllSeeds,main=figTitle,ylab="Heights",xlab='FC windows', labels=F)
		dev.off()
							
		# clustering and cut the tree for windows of each seed
		clustIndxEachSeed <- matrix(nrow=numWinPerSeed,ncol=numSeed,byrow=T)
		starts <- seq(1,totNumWin,by=numWinPerSeed)
		ends <- starts + numWinPerSeed -1
		for (k in seedNum) {
			sprintf('Working on seed %i',k)
			winEachSeed <- winAllSeeds[starts[k]:ends[k],]
			dim(winEachSeed)
			distEachSeed <- dist(winEachSeed, method='euclidean')
			dendroEachSeed <- hclust(distEachSeed,method='ward')
			clustIndxEachSeedTmp <- cutreeDynamic(dendroEachSeed, minClusterSize=2, distM=as.matrix(distEachSeed))
			clustIndxEachSeed[,k] <- clustIndxEachSeedTmp
		
		# plot and save the dendrogram plot
		figName <- sprintf("/Users/zhen.yang/Documents/Zhen_CMI/microstate/fig/lambdaOptimalPerSub/zWinEachSeed_%sCorLasso_%s_%s_10min.png",winType[j], sessionType[i],seed[k])
		png(filename = figName , width = 480, height = 480, units = "px", pointsize = 12,bg = "white")
		figTitle <- sprintf("zWinEachSeed_%sCorLasso_%s_%s_10min",winType[j], sessionType[i],seed[k])
		plot(dendroEachSeed,main=figTitle,ylab="Heights",xlab='FC windows', labels=F)
		dev.off()
			
	}
		# save the index file
		indxFileNameEachSeed <- sprintf("/Users/zhen.yang/Documents/Zhen_CMI/microstate/results/lambdaOptimalPerSub/zWinEachSeed_%sCorLasso_%s_10min.txt", winType[j], sessionType[i])	
 write.table(clustIndxEachSeed,file=indxFileNameEachSeed,row.names=FALSE,col.names=FALSE,qmethod="double")
	}
	}




