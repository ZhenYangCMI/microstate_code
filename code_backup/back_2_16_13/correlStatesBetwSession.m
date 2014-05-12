clear
clc
close all

sessionlLambda=0.12;
session2Lambda=0.11;
resultDir1=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/session1/session1_ward_euclidean/lambda',num2str(session1Lambda)];
resultDir2=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/session2/session2_ward_euclidean/lambda',num2str(session2Lambda)];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/correlTwoSessions/'];

minNumClust=2;
maxNumClust=2;

k=0;

for i=minNumClust:maxNumClust
    k=k+1;
    numClust=i;
    disp(['Working on ', num2str(numClust), 'clusters.'])
    tmp1=load([resultDir1,'/clusterMean_',num2str(numClust),'clusters.mat']);
    tmp2=load([resultDir2,'/clusterMean_',num2str(numClust),'clusters.mat']);
    clustMean1=tmp1.finalMeanWinOfClust;
    clustMeanTransp1=clustMean1';
    clustMean2=tmp2.finalMeanWinOfClust;
    clustMeanTransp2=clustMean2';
    clustMeanTransp=[clustMeanTransp1,clustMeanTransp2];
    corClusters=corrcoef(clustMeanTransp);
    
    figure(k)
    imagesc(corClusters)
    colorbar
    title('Correlations between clusters')
    xlabel('Clusters')
    ylabel('clusters')
    caxis([-1 1])
    saveas(figure(k),[figDir,'CorrelBetwSession_',num2str(numClust),'clusters_ses1Lambda',num2str(session1Lambda),'_ses2Lambda',num2str(session2Lambda),'.png'])
end

% scatterplot the 2 clusters solution

for i=minNumClust:maxNumClust
    numClust=i;
    disp(['Working on ', num2str(numClust), 'clusters.'])
%     tmp1=load([resultDir1,'/clusterMean_',num2str(numClust),'clusters.mat']);
%     tmp2=load([resultDir2,'/clusterMean_',num2str(numClust),'clusters.mat']);
%     clustMean1=tmp1.finalMeanWinOfClust;
%     clustMean2=tmp2.finalMeanWinOfClust;
    figure(1)
    scatter(clustMean1(1,:),clustMean2(2,:))
    xlabel('Session 1 Cluster 1')
    ylabel('Session 2 Cluster 2')
    xlim([-0.06 0.03])
    ylim([-0.14 0.06])
    saveas(figure(1),[figDir,'scatter_ses1State1_ses2State2',num2str(numClust),'clusters_ses1Lambda',num2str(session1Lambda),'_ses2Lambda',num2str(session2Lambda),'.png'])
    figure(2)
    scatter(clustMean1(2,:),clustMean2(1,:))
    xlabel('Session 1 Cluster 2')
    ylabel('Session 2 Cluster 1')
    xlim([-0.06 0.03])
    ylim([-0.14 0.06])
    saveas(figure(2),[figDir,'scatter_ses1State2_ses2State1',num2str(numClust),'clusters_ses1Lambda',num2str(session1Lambda),'_ses2Lambda',num2str(session2Lambda),'.png'])
end



