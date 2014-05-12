clear
clc
close all

analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
figDir=[analyDir, 'fig/645/all_10min/seedTC/']

tmp=load('/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/session1/zWinFullCorLasso_OptimalLambdaPerSub_645_session1.mat');
windows=tmp.zWinFullCorLasso;
winIDList=10:5984:23936
numPlot=length(winIDList);
for i=1:numPlot
    winID=winIDList(i)
    figure(i)
    imagesc(windows(winID,:))
    %colorbar
    caxis([-3 3])
    axis off
    saveas(figure(i),[figDir, 'covVectorDynamic2_seed',num2str(i),'.png'])
end

tmp1=load('/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/session1/z_stationary_FC.mat');
FC=tmp1.z_stationary_FC;
numSeed=4;
subID=10;
for i=1:numSeed
        figure(i)
    imagesc(FC(i,:, subID))
    %colorbar
    caxis([-0.6 0.6])
    axis off
    saveas(figure(i),[figDir, 'covVectorStationary_seed',num2str(i),'sub', num2str(subID),'.png'])
end
