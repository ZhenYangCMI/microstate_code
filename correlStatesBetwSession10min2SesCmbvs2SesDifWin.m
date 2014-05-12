clear
clc
close all

clear
clc

dataLength='all_10min';


resultDir1=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/', dataLength, filesep, 'session1', filesep, 'clustMean'];
resultDir2=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/', dataLength, filesep, '2sessions'];

figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/', dataLength, '/2sessions/'];


    numClust2=4;
    numClust1=5;
    tmp1=load([resultDir1,'/clusterMean_',num2str(numClust1),'clusters_session1_normWin.mat']);
    tmp2=load([resultDir2,'/clusterMean_',num2str(numClust2),'clusters_2sessions_21sub_win136_normWin.mat']);
    tmp3=load([resultDir2,'/clusterMean_',num2str(numClust2),'clusters_2sessions_20sub_win34_normWin.mat']);
    
    clustMean1=tmp1.finalMeanWinOfClust;
    clustMeanTransp1=clustMean1';
    clustMean2=tmp2.finalMeanWinOfClust;
    clustMeanTransp2=clustMean2';
     clustMean3=tmp3.finalMeanWinOfClust;
    clustMeanTransp3=clustMean3';
    
   transp1=[clustMeanTransp1,clustMeanTransp2];
    [corClusters1,pValue1]=corrcoef(transp1);
    
    transp2=[clustMeanTransp1,clustMeanTransp3];
    [corClusters2,pValue2]=corrcoef(transp2);
      
         close all
    figure(1)
    imagesc(corClusters1)
    colorbar
    title('Correlations between clusters of 1session and 2sessions win136')
    xlabel('Clusters (session 1: 1-5; 2sessions 6-9)')
    ylabel('clusters (session 1: 1-5; 2sessions 6-9)')
    caxis([-1 1])
   saveas(figure(1),[figDir,'CorrelBetw2Sesvs1Ses_zWinFullCorLasso_win136.png'])

   figure(2)
    imagesc(corClusters2)
    colorbar
    title('Correlations between clusters of 1session and 2sessions win34')
    xlabel('Clusters (session 1: 1-5; 2sessions 6-9)')
    ylabel('clusters (session 1: 1-5; 2sessions 6-9)')
    caxis([-1 1])
   saveas(figure(2),[figDir,'CorrelBetw2Sesvs1Ses_zWinFullCorLasso_win34.png'])
   
   
   
   
   
   
   
   
   
   
   
   


