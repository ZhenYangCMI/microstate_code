#clear screen
#command+alt+l

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
numSub <- 22
numWinPerSeed <- numWinPerSub * numSub
totNumWin <-numWinPerSub * numSub * numSeed

sessionType <- c('session2')
sessionNum <- 1:length(sessionType)
winType <- c("Full")
winTypeNum <- 1:length(winType)
clustNum <- seq(2,100,by=1)
numClust <- length(clustNum)

for (i in sessionNum) {
	for (j in winTypeNum) {
		baseResults <- sprintf("/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/%s", sessionType[i])
		baseFig <- sprintf("/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/%s", sessionType[i])
		fileName <- file.path(baseResults, sprintf("zWin%sCorLasso_OptimalLambdaPerSub_645_%s.mat", winType[j], sessionType[i]))
		matData <- readMat(fileName)
		# Copy over matrix
		varName <- sprintf("zWin%sCorLasso",winType[j])
		winAllSeeds <- matData[[varName]]
		# Check dimensions
		dim(winAllSeeds)

		# clustering and cut the tree for windows of all seeds 
		distAllSeeds <- dist(winAllSeeds, method='euclidean')
		dendroAllSeeds <- hclust(distAllSeeds,method='ward')
		
		clustIndx=matrix(nrow=totNumWin,ncol=numClust,byrow=F)
		for (m in clustNum) {
		clustIndx[,m-1] <- cutree(dendroAllSeeds, k=clustNum[m-1])
		}
				
		# save the index file
		indxFileName <- file.path(baseResults,sprintf("clustIndx_2to100clusters_zWinAllSeeds_%sCorLasso_%s_10min.txt", winType[j], sessionType[i]))
		write.table(clustIndx,file=indxFileName,row.names=FALSE,col.names=FALSE,qmethod="double")
	}
}
	
	


