# This script will do the hierarchical clustering and using dyamic treecut to estimate the optimal number of clusters
#clear workspace
rm(list = ls(all = TRUE))

# user define the following parameters
covType <- 'GSR'  # covType can be 'compCor' or 'noGSR'
session <- c('session1')
numROI <- c('0100')
winType <- 'winFullCorLasso'
baseResults <- '/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/testing/results/TR645_5min/'
baseFig <- baseResults
winFile <- sprintf("zfeatureFC_%s_%s_%s.mat", winType, numROI, session)
varName <- 'zfeatureWin'

## The followin section dosen't need to change ##


# Load matlab library
library(R.matlab)

# Read in matlab data
fileName <- file.path(baseResults, winFile)
matData <- readMat(fileName)

# Copy over matrix
win <- matData[[varName]]
#win5min <- rbind(win[1:136,],win[273:408,])
# Check dimensions
#dim(win5min)
dim(win)
numWin=dim(win)[1]

# import daynamic tree cut library
library(dynamicTreeCut)

# clustering and cut the tree for windows of all seeds 
dist <- dist(win, method='euclidean')
dendro <- hclust(dist,method='ward')

clustNum <- 1:14
numSolutions <- length(clustNum)

clustIndx <- matrix(nrow=numWin,ncol=numSolutions,byrow=T)
for (i in clustNum) {
			sprintf('Working on %i clusters',i)	
clustIndxEachSolution <- cutree(dendro, k = i)
			clustIndx[,i] <- clustIndxEachSolution
}	
# save the index file
indxFileName <- file.path(baseResults,sprintf("clustIndx_%s_%s_2to14clusters.txt", winType, numROI))
write.table(clustIndx,file=indxFileName,row.names=FALSE,col.names=FALSE,qmethod="double")

# plot and save the dendrogram plot
figName <- file.path(baseFig,sprintf("dendro_%s_%s_.png", winType, numROI))
png(filename = figName , width = 480, height = 480, units = "px", pointsize = 12,bg = "white")
figTitle <- sprintf("denfro_%s_%s", winFile, numROI)
plot(dendro,main=figTitle,ylab="Heights",xlab='feature windows', labels=F)
dev.off()
	

		
							
		

