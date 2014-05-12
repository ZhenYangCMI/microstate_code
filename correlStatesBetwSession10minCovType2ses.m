clear
clc
close all


dataLength='all_10min';

winSize=69;

resultDir1=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/', dataLength, filesep, 'GSR/session1/clustMean/'];
resultDir2=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/', dataLength, filesep, 'GSR/session2/clustMean/'];

resultDir3=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,'/noGSR/session1/'];
resultDir4=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,'/noGSR/session2/'];

resultDir5=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,'/noGSR/2sessions/'];


figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, '/correlDifCovType/'];


numClust1=5;
numClust2=6;
numClust3=4;
numClust4=4;
numClust5=6;

tmp1=load([resultDir1,'clusterMean_',num2str(numClust1),'clusters_session1_normWin.mat']);
tmp2=load([resultDir2,'clusterMean_',num2str(numClust2),'clusters_session2_normWin.mat']);
tmp3=load([resultDir3,sprintf('clusterMean_%dclusters_session1_noGSR_win%d_normWin.mat',numClust3,winSize)]);
tmp4=load([resultDir4,sprintf('clusterMean_%dclusters_session2_noGSR_win%d_normWin.mat',numClust4,winSize)]);
tmp5=load([resultDir5,sprintf('clusterMean_%dclusters_2sessions_noGSR_win%d_normWin.mat',numClust5,winSize)]);

clustMean1=tmp1.finalMeanWinOfClust;
clustMeanTransp1=clustMean1';

clustMean2=tmp2.finalMeanWinOfClust;
clustMeanTransp2=clustMean2';

clustMean3=tmp3.finalMeanWinOfClust;
clustMeanTransp3=clustMean3';

clustMean4=tmp4.finalMeanWinOfClust;
clustMeanTransp4=clustMean4';

clustMean5=tmp5.finalMeanWinOfClust;
clustMeanTransp5=clustMean5';

cmb=[clustMeanTransp1,clustMeanTransp5];
[corClust,pValue]=corrcoef(cmb);


ColorMap=[1,1,0;1,0.9,0;1,0.8,0;1,0.7,0;1,0.6,0;1,0.5,0;0,0.5,1;0,0.6,1;0,0.7,1;0,0.8,1;0,0.9,1;0,1,1;];
ColorMap=flipdim(ColorMap,1);
cmap1 = colorRamp(ColorMap(1:6,:), 32);
cmap2= colorRamp(ColorMap(7:end,:), 32);
ColorMap=vertcat(cmap1,cmap2)

figure(1)
imagesc(corClust)
h=figure(1);
colorbar
caxis([-1 1])
%set(gca,'xTick',0.5:5.5, 'yTick', 0.5:5.5,'XTickLabel',[],'YTickLabel',[]);
set(gca,'LineWidth',2)

%set(figure(1),'Colormap',ColorMap)  % this is used to generate the
%first graph for eiditing the colormap
grid off
set(gca, 'GridLineStyle', '-' )
% xlabel('clusters of GSR: 1-5; clusters of noGSR: 6-10')
% ylabel('clusters of GSR: 1-5; clusters of noGSR: 6-10')
% title('Correlations between clusters of GSR and noGSR')
% set(h,'Position',[1,1,800,665])
saveas(figure(1),[figDir,'CorrelGSRSes1AndnoGSR2sessions_zWinFullCorLasso.png'])







