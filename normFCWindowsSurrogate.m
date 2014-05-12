clear
clc
close all

sessionList={'session1'};
numSession=length(sessionList)
numSurrogate=100;

for j=1:numSurrogate
    for i=1:numSession
        session=char(sessionList{i});
        resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/',session,'/surrogate/'];
        outputDir=[resultDir,'zNormWindows/'];
        tmp1 =load([resultDir, 'unNormWindows/winFullCorLasso_OptimalLambdaPerSub_645_',session,'_surrogate_',num2str(j),'.mat']);
        winFullCorLasso=tmp1.winFullCorLasso;
        tmp2=winFullCorLasso';
        numROI=size(tmp2,1);
        tmp3=(tmp2-repmat(mean(tmp2),numROI,1))./(repmat(std(tmp2),numROI,1));
        zWinFullCorLasso=tmp3';
        save([outputDir,'zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'_surrogate_',num2str(j),'.mat'],'zWinFullCorLasso')
    end
end