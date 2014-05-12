# This script will do the hierarchical clustering and using dyamic treecut to estimate the optimal number of clusters
#clear workspace
rm(list = ls(all = TRUE))

# user define the following parameters
covType <- 'GSR'  # covType can be 'compCor' or 'noGSR'
session <- c('session1')
numWinPerSub <- 272
totNumWin <-5984
numSub <- 22
winType <- 'winFullCorLasso'

baseResults <- sprintf("/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/results/%s/%s", covType, session)
baseOutput <- sprintf("/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/results/%s/%s/eachSub", covType, session)
baseFig <- sprintf("/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/figs/%s/clustering/%s", covType, session)
winFile <- sprintf("zfeatureFC_%s_%s.mat", winType, session)
varName <- 'zfeatureWin'

## The followin section dosen't need to change ##


# Load matlab library
library(R.matlab)

# Read in matlab data
fileName <- file.path(baseResults, winFile)
matData <- readMat(fileName)
#subList <- readMat("/home/data/Projects/microstate/NKITRT_SubID.mat")

# Copy over matrix
win <- matData[[varName]]


# Check dimensions
dim(win)
dim(subList)

sub <- 1:numSub

# import daynamic tree cut library
library(dynamicTreeCut)

clustIndxEachSub <- matrix(nrow=numWinPerSub,ncol=numSub,byrow=T)
		starts <- seq(1,totNumWin,by=numWinPerSub)
		ends <- starts + numWinPerSub -1
		for (k in sub) {
			sprintf('Working on sub %i',k)
			winEachSub <- win[starts[k]:ends[k],]
			dim(winEachSub)
			distEachSub <- dist(winEachSub, method='euclidean')
			dendroEachSub <- hclust(distEachSub,method='ward')
			clustIndxEachSubTmp <- cutreeDynamic(dendroEachSub, minClusterSize=2, distM=as.matrix(distEachSub))
			clustIndxEachSub[,k] <- clustIndxEachSubTmp
		
		# plot and save the dendrogram plot
		figName <- file.path(baseFig,sprintf("zFeatureWinEachSub_Fullcor_%i_%s.png",sub[k],session))
		png(filename = figName , width = 480, height = 480, units = "px", pointsize = 12,bg = "white")
		figTitle <- sprintf("zFeatureWinEachSub_Fullcor_%i_%s.png",sub[k],session)
		plot(dendroEachSub,main=figTitle,ylab="Heights",xlab='feature windows', labels=F)
		dev.off()
		}
		# save the index file
		indxFileNameEachSub <- file.path(baseOutput,sprintf("clustIndx_zFeatureWinEachSub_%s_%s.txt", winType, session))	
 write.table(clustIndxEachSub,file=indxFileNameEachSub,row.names=FALSE,col.names=FALSE,qmethod="double")
	

		

		
		
							
		

