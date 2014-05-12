% This script normalize the iFC windows before clustering

clear
clc
close all


session='2sessions';
covType='noGSR'; % 'compCor', 'GSR', or 'noGSR'
winSize=69;

resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/',covType, '/', session,'/'];
tmp1 =load([resultDir, 'winFullCorLasso_OptimalLambdaPerSub_645_session1win', num2str(winSize), '_', covType, '.mat']);
winFullCorLasso1=tmp1.winFullCorLasso;
tmp2 =load([resultDir, 'winFullCorLasso_OptimalLambdaPerSub_645_session2win', num2str(winSize), '_', covType, '.mat']);
winFullCorLasso2=tmp2.winFullCorLasso;

winFullCorLasso=vertcat(winFullCorLasso1, winFullCorLasso2);

tmp=winFullCorLasso';
numROI=size(tmp,1);
tmp3=(tmp-repmat(mean(tmp),numROI,1))./(repmat(std(tmp),numROI,1));
zWinFullCorLasso=tmp3';
save([resultDir, 'zwinFullCorLasso_OptimalLambdaPerSub_645_',session,'win', num2str(winSize), '_', covType, '.mat'],'zWinFullCorLasso')


