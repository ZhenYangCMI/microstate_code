clear
clc

numSub=22;
covType='GSR';
session='session2';
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
resultDir=[analyDir,'fullBrainFC/results/',covType, '/', session,'/test/'];
figDir=[analyDir,'fullBrainFC/figs/', covType, '/'];

tmp1=load([resultDir,'featureFC_winFullCor_',session,'.mat']);
tmp2=load([resultDir,'featureFC_winFullCorLasso_',session,'.mat'])
tmp3=load([resultDir,'featureFC_winPartialCorLasso_',session,'.mat'])
winFullCor=tmp1.featureWin;
winFullCorLasso=tmp2.featureWin;
winPartialCorLasso=tmp3.featureWin;

numWinPerSub=size(winFullCor, 1)/numSub;
% sub are orderd according to low, mid, and high motion

subList=[20 17 18];


for j=1:length(subList)
    close all
    sub=subList(j)
    concate=[mean(winFullCor(numWinPerSub*(j-1)+1:numWinPerSub*j, :))', mean(winFullCorLasso(numWinPerSub*(j-1)+1:numWinPerSub*j, :))', mean(winPartialCorLasso(numWinPerSub*(j-1)+1:numWinPerSub*j, :))'];
    [r, p]=corrcoef(concate)
    
    winType={'winFullCor', 'winFullCorLasso', 'winPartialCorLasso'};
    
    t=0
    for i=1:length(winType)
        x=char(winType{i})
        for j=1:length(winType)
            y=char(winType{j})
            t=t+1
            figure(1)
            subplot(3,3, t)
            scatter(concate(:,i), concate(:, j))
            xlim([-4 4])
            ylim([-4 4])
            lsline
        end
    end
    saveas(figure(1), [figDir, 'scatter_', x, 'vs_', y, 'sub', num2str(sub), '_', session, '.png'])
end