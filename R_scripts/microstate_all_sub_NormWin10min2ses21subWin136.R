#clear screen
ctl+alt+l

#clear workspace
rm(list = ls(all = TRUE))

###
# Load Data
###

# Load matlab library
library(R.matlab)

# Read in matlab data

winSize <-136
sessionType <- c('2sessions')
sessionNum <- 1:length(sessionType)
winType <- c("Full")
winTypeNum <- 1:length(winType)

seed <- c('seed1','seed2','seed3','seed4' )
seedNum <- 1:length(seed)
i <- 1
j <- 1

for (i in sessionNum) {
	for (j in winTypeNum) {

		baseResults <- sprintf("/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/%s", sessionType[i])
		baseFig <- sprintf("/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/%s", sessionType[i])
		fileName <- file.path(baseResults, sprintf("zWinFullCorLasso_OptimalLambdaPerSub_645_2sessions_21sub_win%d.mat", winSize))
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
	indxFileName <- file.path(baseResults,sprintf("clustIndxNormWinAllSeeds_%sCorLasso_%s_10min_21sub_win%d.txt", winType[j], sessionType[i], winSize))
	write.table(clustIndx,file=indxFileName,row.names=FALSE,col.names=FALSE,qmethod="double")
	
	# plot and save the dendrogram plot
		figName <- file.path(baseFig,sprintf("zWinAllSeeds_%sCorLasso_%s_10min_21sub_win%d.png",winType[j], sessionType[i], winSize))
		png(filename = figName , width = 480, height = 480, units = "px", pointsize = 12,bg = "white")
		figTitle <- sprintf("zWinAllSeeds_%sCorLasso_%s_10min_21sub_win%d",winType[j], sessionType[i], winSize)
		plot(dendroAllSeeds,main=figTitle,ylab="Heights",xlab='FC windows', labels=F)
		dev.off()
		}
		}
}
		
		
							
		# clustering and cut the tree for windows of each seed
numWinPerSub <-250 # 250 for win34 (22s); 284 for win136 (88s)
numSeed <- 4
numSub <- 21
numWinPerSeed <- numWinPerSub * numSub
totNumWin <-numWinPerSub * numSub * numSeed
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
		figName <- file.path(baseFig,sprintf("zWinEachSeed_%sCorLasso_%s_%s_10min.png",winType[j], sessionType[i],seed[k]))
		png(filename = figName , width = 480, height = 480, units = "px", pointsize = 12,bg = "white")
		figTitle <- sprintf("zWinEachSeed_%sCorLasso_%s_%s_10min",winType[j], sessionType[i],seed[k])
		plot(dendroEachSeed,main=figTitle,ylab="Heights",xlab='FC windows', labels=F)
		dev.off()
		}
		# save the index file
		indxFileNameEachSeed <- file.path(baseResults,sprintf("zWinEachSeed_%sCorLasso_%s_10min.txt", winType[j], sessionType[i]))	
 write.table(clustIndxEachSeed,file=indxFileNameEachSeed,row.names=FALSE,col.names=FALSE,qmethod="double")
	}
	}




