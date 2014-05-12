function [indx,numclust,clustTree,cophenetCoef,inconsistCoef]=hierarchicalCluster...
    (corType, cutoffThresh, winAllSubAllSeedFullCor, winAllSubAllSeedPartialCor)
% This function returns the clustering indx for each window
% and the centroids for each cluster
%input
%1. sessionList: cell, e.g. {'session1','session2'}
%2. corType: string, 'full' or 'partial', one type at a time
%3. cutoffThresh: numerical, e.g. 0.1


if strcmp('full',corType)
    winForClustering=winAllSubAllSeedFullCor;
elseif strcmp('partial',corType)
    winForClustering=winAllSubAllSeedPartialCor;
else
    disp('Incorrect input for corType.')
end
distBtwWin= pdist(winForClustering,'correlation');
clustTree = linkage(distBtwWin,'ward');
cophenetCoef = cophenet(clustTree,distBtwWin)
indx = cluster(clustTree,'cutoff',cutoffThresh,'criterion','distance' );
numclust=max(indx)
inconsistCoef=inconsistent(clustTree);
dendrogram(clustTree, 0)


