clear
clc
close all

clear
clc

dataLength='all_10min';

resultDir1=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/', dataLength, '/GSR/session1/clustMean/'];
resultDir2=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/', dataLength, '/GSR/session2/clustMean/'];
%figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/', dataLength, '/correlTwoSessions/'];


    
    numClust1=5;
    numClust2=6;
    tmp1=load([resultDir1,'/clusterMean_',num2str(numClust1),'clusters_session1_normWin.mat']);
    tmp2=load([resultDir2,'/clusterMean_',num2str(numClust2),'clusters_session2_normWin.mat']);
    clustMean1=tmp1.finalMeanWinOfClust;
    clustMeanTransp1=clustMean1';
    clustMean2=tmp2.finalMeanWinOfClust;
    clustMeanTransp2=clustMean2';
    clustMeanTransp=[clustMeanTransp1,clustMeanTransp2];
    [corClusters,pValue]=corrcoef(clustMeanTransp);
    
    figure(1)
    imagesc(corClusters)
    colorbar
    title('Correlations between clusters')
    xlabel('Clusters')
    ylabel('clusters')
    caxis([-1 1])
    saveas(figure(1),[figDir,'CorrelBetwTwoSessions_zWinFullCorLasso.png'])



