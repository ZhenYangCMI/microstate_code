#clear workspace
rm(list = ls(all = TRUE))

#clear screen


# Load Data
###

# Load matlab library
library(R.matlab)

# Read in matlab data
base <- "/home/data/Projects/microstate/DPARSF_preprocessed/results/645/session1"
matdata <- readMat(file.path(base, "winAllSubAllSeed_partialCor_645_session1_0.12.mat"))

# Copy over matrix
mat <- matdata$winAllSubAllSeedPartialCor
rm(matdata)

# Check dimensions
dim(mat)

# See sample data
mat[1:5,1:5]


###
# Hierarchical Clustering
###

# Calculate euclidean distance
dmat <- dist(mat)

# Compute dendrogram / hierarchical clustering
dendro <- hclust(dmat, "ward")
write(dendro,file=(file.path(base,"hierachical clustering R"))

# Cluster
library(dynamicTreeCut)
clust <- cutreeDynamic(dendro, minClusterSize=2, distM=as.matrix(dmat))

# If you want to plot the dendrogram
plot(dendro)

# See breakdown of clusters
table(clust)

# Can do simpler clustering
members <- cutree(dendro, k=10) # gives you 10 clusters


