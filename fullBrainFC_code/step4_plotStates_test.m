%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script do the following
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

% load subList
subList={'3808535','8574662'};
numSub=length(subList);
numROI='0050';
session='session1'
clustMethod='kmeans'

resultDir='/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/testing/results/TR2500_5min/'
figDir=[resultDir, '/figs/'];
% load the index file
if strcmp(clustMethod, 'hierachical')
clustIndx=load([resultDir, 'clustIndx_winFullCorLasso_', numROI, '_2to14clusters.txt']);
else
clustIndx=load([resultDir, 'clustIndx_', numROI, '_1to20200rep.mat']);
clustIndx=clustIndx.index;
end

win=load([resultDir, 'zfeatureFC_winFullCorLasso_', numROI, '_', session, '.mat'])
winData=win.zfeatureWin;
numFeature=size(winData, 2);

maxNumClust=14

for i=1:maxNumClust
    numClust=i
    % plot the cluster mean
    figure(i)
    clustMean=zeros(numClust, numFeature);
    for j=1:numClust
        clust=winData(find(clustIndx(:,i)==j), :);
        clustMean(j,:)=mean(clust);
        clustFC=squareform(clustMean(j,:));
        if i> 7
            subplot(7,2,j)
        else
            subplot(7,1,j)
        end
        imagesc(clustFC)
        colorbar
        caxis([-1 1])
        %axis square
    end
    saveas(figure(i), [figDir, session, '_', numROI, '_stateMap_', num2str(i), 'clusters_', clustMethod, '.png'])
end







