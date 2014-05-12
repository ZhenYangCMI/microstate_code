clear
clc
close all

clear
clc

dataLength='all_10min';

resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/', dataLength, filesep, '2sessions'];

figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/', dataLength, '/2sessions/correlOddEven/'];


    
    numClust=3;
    
    tmp1=load([resultDir,'/clusterMean_',num2str(numClust),'clusters_2sessions_odd_normWin.mat']);
    tmp2=load([resultDir,'/clusterMean_',num2str(numClust),'clusters_2sessions_even_normWin.mat']);
    clustMean1=tmp1.finalMeanWinOfClust;
    clustMeanTransp1=clustMean1';
    clustMean2=tmp2.finalMeanWinOfClust;
    clustMeanTransp2=clustMean2';
    clustMeanTransp=[clustMeanTransp1,clustMeanTransp2];
    [corClusters,pValue]=corrcoef(clustMeanTransp);
    
    figure(1)
    imagesc(corClusters)
    colorbar
    title('Correlations between clusters of even and odd windows')
    xlabel('Clusters (odd: 1-4; even 5-8)')
    ylabel('clusters (odd: 1-4; even 5-8)')
    caxis([-1 1])
    saveas(figure(1),[figDir,'CorrelBetwOddEven_zWinFullCorLasso.png'])



