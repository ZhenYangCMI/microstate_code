clear
clc
close all

clear
clc

dataLength='all_10min';
session='session1';

resultDir1=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength, filesep, session,'/clustMean'];
resultDir2=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/', dataLength, filesep, session, '/surrogate/clustMean'];

figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/', dataLength,filesep, session,'/surrogateData/'];



numClust=5;
numSurrogate=100;

% ColorMap=[1,1,0;1,0.9,0;1,0.8,0;1,0.7,0;1,0.6,0;1,0.5,0;0,0.5,1;0,0.6,1;0,0.7,1;0,0.8,1;0,0.9,1;0,1,1;];
% ColorMap=flipdim(ColorMap,1);
% cmap1 = colorRamp(ColorMap(1:6,:), 32);
% cmap2= colorRamp(ColorMap(7:end,:), 32);
% ColorMap=vertcat(cmap1,cmap2)

 %load('MyColormapsCorrel','mycmap')

finalCorClusters=zeros(2*numClust, 2*numClust,numSurrogate);

for i=1:numSurrogate
    disp(['Working on surrogate ', num2str(i)])
    
    tmp1=load([resultDir1,'/clusterMean_',num2str(numClust),'clusters_', session, '_normWin.mat']);
    tmp2=load([resultDir2,'/clusterMean_',num2str(numClust),'clusters_session1_normWin_surrogate_', num2str(i),'.mat']);
    clustMean1=tmp1.finalMeanWinOfClust;
    clustMeanTransp1=clustMean1';
    clustMean2=tmp2.finalMeanWinOfClust;
    clustMeanTransp2=clustMean2';
    clustMeanTransp=[clustMeanTransp1,clustMeanTransp2];
    [corClusters,pValue]=corrcoef(clustMeanTransp);
    finalCorClusters(:,:,i)=corClusters;
    
end

tmp=reshape(finalCorClusters,[],numSurrogate);
tmp=tmp';
meanCor=mean(tmp);
finalMean=reshape(meanCor, 2*numClust,2*numClust);

figure(1)
imagesc(finalMean)
colorbar
title('Mean correlations between clusters of True and surrogated data')
xlabel('Clusters (True: 1-5; Surrogate: 6-10)')
ylabel('clusters (True: 1-5; Surrogate: 6-10)')
caxis([-1 1])
%set(figure(1),'Colormap',ColorMap)
saveas(figure(1),[figDir,'CorrelBetwTrueVsSurrogate_session1_zWinFullCorLasso.png'])

figure(2)
data1=finalCorClusters(6:end,1:5,:);
dataPlot1=reshape(data1,1,[]);
hist(dataPlot1, -1:0.1:1)
title('Correl: True and Surrogate')
ylim([0,500])
xlim([-1,1])
saveas(figure(2),[figDir,'session1_corDistribution_TrueAndSurrogate.png'])

figure(3)
data2=finalCorClusters(6:end,6:end,:);
dataPlot2=reshape(data2,1,[]);
hist(dataPlot2, -1:0.1:1)
title('Correl: Surrogate and Surrogate')
ylim([0 500])
xlim([-1 1])
saveas(figure(3),[figDir,'session1_corDistribution_SurrogateAndSurrogate.png'])



