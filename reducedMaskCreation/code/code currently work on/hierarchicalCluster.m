function [indx,clustTree,cophenetCoef, inconsistCoef]=hierarchicalCluster(winForClustering, cutoffThresh)
% This function returns the clustering indx for each window 
% and the centroids for each cluster

distBtwWin = pdist(winForClustering,'correlation'); 
clustTree = linkage(distBtwWin,'single'); 
cophenetCoef = cophenet(clustTree,distBtwWin);
indx = cluster(clustTree,'cutoff',cutoffThresh,'criterion','distance' )
inconsistCoef=inconsistent(clustTree);