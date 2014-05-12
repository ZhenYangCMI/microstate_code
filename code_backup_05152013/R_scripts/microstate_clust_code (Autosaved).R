#clear workspace
rm(list = ls(all = TRUE))

#clear screen
command+alt+l

###
# Load Data
###

# Load matlab library
library(R.matlab)

# Read in matlab data
base <- "/Users/zhen.yang/Desktop"
base <- "/home/data/Projects/microstate/DPARSF_preprocessed/results/645/session1"
matdata <- readMat(file.path(base, "winAllSubAllSeed_partialCor_645_session1_0.12.mat"))

# Copy over matrix
mat <- matdata$winAllSubAllSeedPartialCor
rm(matdata)

# Check dimensions
dim(mat)

library(dynamicTreeCut)

clusts <- matrix(nrow=2816,ncol=4,byrow=T)

xis <- 1:4
starts <- seq(1,11264,by=2816)
ends <- starts + 2816 - 1
for (i in xis) {
	print(i)
	seed <- mat[starts[i]:ends[i],]
	dim(seed)
	dData <- dist(seed, method='euclidean')
	dendroH <- hclust(dData,method='ward')
	cluststemp <- cutreeDynamic(dendroH, minClusterSize=2, distM=as.matrix(dData))
	clusts[,i] <- cluststemp
}

clusts[[1]]

clusts[[2]]
clusts[[3]]


# See sample data
seed1 <- mat[1:2816,]


###
# Hierarchical Clustering
###
dData <- dist(mat, method='euclidean')
# Calculate euclidean distance and apply the agglorative nesting (hierachical clustering)
# the dist and hclust function is different from matlab linkage function, hclust used a lance-williams dissimilarity update formular, also the distance computed is . use agnes, if the method is not "flexible", the the lance-williams formula will not be applied
dendro <-agnes(x = seed1, diss = F, method = "ward",metric='euclidean')
dendroH <- hclust(dData,method='ward')

# Cluster
library(dynamicTreeCut)
clust <- cutreeDynamic(dendroH, minClusterSize=2, distM=as.matrix(dData))

# If you want to plot the dendrogram
plot(dendro)

# See breakdown of clusters
table(clust)

# Can do simpler clustering
members <- cutree(dendro, k=10) # gives you 10 clusters


