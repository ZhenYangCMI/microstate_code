clear
clc
close all

clear
clc

dataLength='all_10min';

resultDir1=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/', dataLength, filesep, '2sessions'];
resultDir2=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/', dataLength, filesep, 'session1', filesep, 'clustMean'];

figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/', dataLength, '/2sessions/'];


    
    numClust=5;
    
    tmp1=load([resultDir1,'/clusterMean_',num2str(numClust),'clusters_2sessions_21sub_normWin.mat']);
    tmp2=load([resultDir2,'/clusterMean_',num2str(numClust),'clusters_session1_normWin.mat']);
    clustMean1=tmp1.finalMeanWinOfClust;
    clustMeanTransp1=clustMean1';
    clustMean2=tmp2.finalMeanWinOfClust;
    clustMeanTransp2=clustMean2';
    clustMeanTransp=[clustMeanTransp1,clustMeanTransp2];
    [corClusters,pValue]=corrcoef(clustMeanTransp);
    
    [pID,pN] = FDR(pValue,0.05)
    
    figure(1)
    imagesc(corClusters)
    colorbar
    title('Correlations between clusters of 2sessions and 1session')
    xlabel('Clusters (2sessions: 1-5; 1session 6-10)')
    ylabel('clusters (2sessions: 1-5; 1session 6-10)')
    caxis([-1 1])
   % saveas(figure(1),[figDir,'CorrelBetw2Sesvs1Ses_zWinFullCorLasso.png'])



