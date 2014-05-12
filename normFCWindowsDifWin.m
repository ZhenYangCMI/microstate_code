clear
clc
close all

% sessionList={'session1','session2'};
sessionList={'session2'};
winSize=34; % test 34 or 136 TRs corresponding to 22 and 88 sec
numSession=length(sessionList)
for i=1:numSession
    session=char(sessionList{i});
    resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/',session,'/'];
    tmp1 =load([resultDir, 'winFullCorLasso_OptimalLambdaPerSub_645_',session,'win',num2str(winSize),'.mat']);
    winFullCorLasso=tmp1.winFullCorLasso;
    tmp2=winFullCorLasso';
    numROI=size(tmp2,1);
    tmp3=(tmp2-repmat(mean(tmp2),numROI,1))./(repmat(std(tmp2),numROI,1));
    zWinFullCorLasso=tmp3';
    save([resultDir,'zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'win',num2str(winSize),'.mat'],'zWinFullCorLasso')
end