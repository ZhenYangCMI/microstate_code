% plot the num of cluster per sub
clear
clc
close all
session='session1'
numSub=22
win='FullcorLasso'
numClustPerSub=zeros(numSub,1);
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/figs/GSR/']
clustIndxEachSub=load(['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/results/GSR/', session, '/eachSub/clustIndx_zFeatureWinEachSub_', win, '.txt']);
for i=1:numSub
    tmp=length(unique(clustIndxEachSub(:, i)));
    numClustPerSub(i,1)=tmp;
end
figure(1)
plot(numClustPerSub)
ylim([12 22])
saveas(figure(1), [figDir, session, '_', win, '_numClustPerSub.png'])
