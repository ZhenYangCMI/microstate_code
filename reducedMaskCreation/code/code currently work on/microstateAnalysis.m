% This script will run the microstate analysis

clear
clc
close all

%TRList={'645','2500'};
%SesList={'session1','session2'};
%load(***)
%subList=SubID;

TRList={'645'};
sessionList={'session1'};
subList={'0021002', '0021018'};

numSub=length(subList);
numSession=length(sessionList);
numTR=length(TRList);
numSeed=4;
numROI=156;

% corType=1(full correlation) or 2 (regularized inverse covariance) 
corType=1
[winAllSubAllSeed]=FCWinCreation(TRList,sessionList,subList, corType);

% Estimate cluster number
winForClustering=winAllSubAllSeed;
cutoffThresh=0.1;
[indx,clustTree,numCluster,cophenetCoef,inconsistCoef]=hierarchicalCluster(winForClustering,cutoffThresh);
figure(1)
dendrogram(clustTree)
saveas(figure(1),[figDir,'/dendrogram.png'])

%final clustering
indxFinal = cluster(clustTree,'maxclust',2)

% percentage of windows per cluster




