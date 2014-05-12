%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script do the following
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

% load subList
subList={'3808535','8574662'};
numSub=length(subList);
numROIList={'0050', '0100', '0200'};
session='session1'
TRList={'645', '2500'};

for m=1:length(numROIList)
    numROI=char(numROIList{m})
    for n=1:length(TRList)
        TR=char(TRList{n})
        close all
        resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/testing/results/TR', TR, '_5min/']
        figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/testing/figs/'];
        % load the index file
        
        clustIndx1=load([resultDir, 'clustIndx_winFullCorLasso_', numROI, '_2to14clusters.txt']);
        
        clustIndx2=load([resultDir, 'clustIndx_1to20clusters_500replications_', numROI, 'ROIs.mat']);
        clustIndx2=clustIndx2.index;
        
        
        win=load([resultDir, 'zfeatureFC_winFullCorLasso_', numROI, '_', session, '.mat'])
        winData=win.zfeatureWin;
        numFeature=size(winData, 2);
        
        % maxNumClust=14
        % t=0;
        % figure(1)
        % for i=3:maxNumClust
        %     numClust=i
        %     t=t+1;
        %     % plot the cluster mean
        %
        %     clustMedian1=zeros(numClust, numFeature);
        %     clustMedian2=zeros(numClust, numFeature);
        %     for j=1:numClust
        %         clust1=winData(find(clustIndx1(:,i)==j), :);
        %         clustMedian1(j,:)=median(clust1);
        %
        %         clust2=winData(find(clustIndx2(:,i)==j), :);
        %         clustMedian2(j,:)=median(clust2);
        %     end
        %     clustMedian1=clustMedian1';
        %     clustMedian2=clustMedian2';
        %
        %     RHO=corr(clustMedian1, clustMedian2);
        %     subplot(3,4,t)
        %     imagesc(RHO)
        %     colorbar
        %     caxis([-1 1])
        %     title([num2str(i), 'clusters'])
        % end
        % saveas(figure(1), [figDir, 'clustMean_betweenClustMethods_', numROI, '_2 to', num2str(maxNumClust), 'clusters_TR', TR, '.png'])
        
        
        maxNumClust=14
        t=0;
        figure(2)
        for i=3:maxNumClust
            numClust=i
            t=t+1;
            % plot the cluster mean
            
            clustMedian1=zeros(numClust, numFeature);
            
            for j=1:numClust
                clust1=winData(find(clustIndx1(:,i)==j), :);
                clustMedian1(j,:)=median(clust1);
            end
            clustMedian1=clustMedian1';
            
            RHO=corr(clustMedian1);
            subplot(3,4,t)
            imagesc(RHO)
            colorbar
            caxis([-1 1])
            title([num2str(i), 'clusters'])
        end
        saveas(figure(2), [figDir, 'clustMean_hierachical_', numROI, '_2 to', num2str(maxNumClust), 'clusters_TR', TR, '.png'])
        
        maxNumClust=14
        t=0;
        figure(3)
        for i=3:maxNumClust
            numClust=i
            t=t+1;
            % plot the cluster mean
            
            clustMedian2=zeros(numClust, numFeature);
            
            for j=1:numClust
                clust2=winData(find(clustIndx2(:,i)==j), :);
                clustMedian2(j,:)=median(clust2);
            end
            clustMedian2=clustMedian2';
            
            RHO=corr(clustMedian2);
            subplot(3,4,t)
            imagesc(RHO)
            colorbar
            caxis([-1 1])
            title([num2str(i), 'clusters'])
        end
        saveas(figure(3), [figDir, 'clustMean_kmeans_', numROI, '_2 to', num2str(maxNumClust), 'clusters_TR', TR, '.png'])
        
    end
end

