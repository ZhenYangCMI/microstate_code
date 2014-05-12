clear
clc
close all

numWinPerSeed=5984;
numSeed=4;
figDir='/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/GSR/session1/';
dataDir='/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/GSR/session1/';

% load the iFC windows of all seeds
a=load([dataDir, 'zWinFullCorLasso_OptimalLambdaPerSub_645_session1.mat'])
tmp=a. zWinFullCorLasso;

% load the index file
indx=load([dataDir, 'clustIndxNormWinAllSeeds_FullCorLasso_session1_10min.txt']);

% recode state 1 and state 2 to match the figs in the paper
indxRecode=indx;
indxRecode(indx==1)=2;
indxRecode(indx==2)=1;
numClust=length(unique(indxRecode));
numROI=size(tmp, 2);
numWin=size(tmp, 1);

clustMean=zeros(numClust, numROI);
for i=1:numClust
winState=tmp(find(indxRecode==i), :);
avg=mean(winState);
clustMean(i, :)=avg;
end

corAllWin=zeros(numWin, 1);
for j=1:numWin
cor=corr(tmp(j, :)', clustMean(indxRecode(j), :)');
corAllWin(j,1)=cor;
end


y=zeros(numWinPerSeed*1, 1);
for i=1:numSeed
figure(i)
stateIndx=indxRecode(1+numWinPerSeed*(i-1):numWinPerSeed*i);
plot(stateIndx, '-*')
ylim([-1 6])
seedCor=corAllWin(1+numWinPerSeed*(i-1):numWinPerSeed*i);
hold on
plot(seedCor, 'r')
hold on
plot(y, '--k')
saveas(figure(i), [figDir, 'TCstateTransitCorWithStateMean_seed', num2str(i), '.png'])
end

close all
subID=21;
figure
for i=1:numSeed
    seed=i;
subplot(2,2,i)
plot(indxRecode(1+272*(subID-1)+5984*(i-1):272*subID+5984*(i-1)))
hold on
plot(corAllWin(1+272*(subID-1)+5984*(i-1):272*subID+5984*(i-1)), 'r')
ylim([-1 6])
title(['seed', num2str(i)])
end

