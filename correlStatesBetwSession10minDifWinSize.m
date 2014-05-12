clear
clc
close all


dataLength='all_10min';
session='session1';

resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/', dataLength, filesep, 'session1/clustMean/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep, session,'/corrWithinSession/'];

winSize1=34;
winSize2=136;

numClust=5;
numClust1=6;
numClust2=5;
tmp=load([resultDir,'clusterMean_',num2str(numClust),'clusters_session1_normWin.mat']);
tmp1=load([resultDir,'clusterMean_',num2str(numClust1),'clusters_session1_normWin_win', num2str(winSize1),'.mat']);
tmp2=load([resultDir,'clusterMean_',num2str(numClust2),'clusters_session1_normWin_win', num2str(winSize2),'.mat']);
clustMean=tmp.finalMeanWinOfClust;
clustMeanTransp=clustMean';

clustMean1=tmp1.finalMeanWinOfClust;
clustMeanTransp1=clustMean1';

clustMean2=tmp2.finalMeanWinOfClust;
clustMeanTransp2=clustMean2';


cmb1=[clustMeanTransp,clustMeanTransp1];
cmb2=[clustMeanTransp,clustMeanTransp2];
[corClust1,pValue1]=corrcoef(cmb1);
[corClust2,pValue2]=corrcoef(cmb2);


ColorMap=[1,1,0;1,0.9,0;1,0.8,0;1,0.7,0;1,0.6,0;1,0.5,0;0,0.5,1;0,0.6,1;0,0.7,1;0,0.8,1;0,0.9,1;0,1,1;];
ColorMap=flipdim(ColorMap,1);
cmap1 = colorRamp(ColorMap(1:6,:), 32);
cmap2= colorRamp(ColorMap(7:end,:), 32);
ColorMap=vertcat(cmap1,cmap2)

figure(1)
imagesc(corClust1)
h=figure(1);
colorbar
caxis([-1 1])
%set(gca,'xTick',0.5:5.5, 'yTick', 0.5:5.5,'XTickLabel',[],'YTickLabel',[]);
set(gca,'LineWidth',2)

%set(figure(1),'Colormap',ColorMap)  % this is used to generate the
%first graph for eiditing the colormap
grid off
set(gca, 'GridLineStyle', '-' )
xlabel('clusters of win44s: 1-5; clusters of win22: 6-11')
ylabel('clusters of win44s: 1-5; clusters of win22: 6-11')
title('Correlations between clusters of win44s and win22s')
set(h,'Position',[1,1,800,665])
saveas(figure(1),[figDir,'Correlwin69AndWin34_zWinFullCorLasso.png'])

figure(2)
imagesc(corClust2)
h=figure(2);
colorbar
caxis([-1 1])
%set(gca,'xTick',0.5:5.5, 'yTick', 0.5:5.5,'XTickLabel',[],'YTickLabel',[]);
set(gca,'LineWidth',2)

grid off
set(gca, 'GridLineStyle', '-' )
xlabel('clusters of win44s: 1-5; clusters of win88: 6-11')
ylabel('clusters of win44s: 1-5; clusters of win88: 6-11')
title('Correlations between clusters of win44s and win88s')
set(h,'Position',[1,1,800,665])
saveas(figure(2),[figDir,'Correlwin69AndWin136_zWinFullCorLasso.png'])






