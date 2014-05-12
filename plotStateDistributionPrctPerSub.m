clear
clc
close all

numWinPerSeed=5984;
numSeed=4;
figDir='/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/GSR/session1/';
dataDir='/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/GSR/session1/';

% load the index file
indx=load([dataDir, 'clustIndxNormWinAllSeeds_FullCorLasso_session1_10min.txt']);

% recode state 1 and state 2 to match the figs in the paper
indxRecode=indx;
indxRecode(indx==1)=2;
indxRecode(indx==2)=1;
numClust=length(unique(indxRecode));

prctWin=zeros(numClust,1);
for i=1:numClust
numWinPerState=length(find(indxRecode==i));
prctWin(i,1)=numWinPerState/(numWinPerSeed*numSeed)*100;
end

lineWidth=2
figure(1)
h=bar(prctWin, 0.45) % TCPrctOverlapStateSum used 0.6, other 2 measures used 0.4 to make them the same width
    set(gca,'LineWidth',lineWidth)
        set(h,'LineWidth',lineWidth,'Facecolor',[0 0 1],'EdgeColor',[0 0 1])
   set(gca,'xTick',1:5, 'yTick', 0:5:30, 'box','off', 'xTickLabel',[],'yTickLabel',[]);
        ylim([0 30])
        saveas(figure(1), [figDir, 'prctWinPerState.png'])

[seedPrct, subPrct, subSeedPrct] = DoQuickStats(indxRecode,numClust)

% the colormap used are manually generated using inversed hot color, then change the black to dark red.
dataPlot=subPrct'*100;
figure(2)
imagesc(dataPlot)
caxis([0 60])
hcb=colorbar
set(hcb,'YTick',[0, 20, 40, 60])
colormap(a.mycmap)
set(gca, 'YTick', [], 'XTick', [])
saveas(figure(2), [figDir, 'prctTimePerSubInEachState.png'])
