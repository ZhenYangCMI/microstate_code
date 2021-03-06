
clear all
clc

covType='noGSR'
session='2sessions'
winSize=69;
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/',covType, '/', session,'/'];
a=load([resultDir, 'zwinFullCorLasso_OptimalLambdaPerSub_645_',session,'win', num2str(winSize), '_', covType, '.mat'])
b=a.zWinFullCorLasso;

% x=randsample(22,1)
% Get an answer of 18, subject 18 will be excluded from this analysis due
% to R object limit
% subID=x;
subID=18;
numROI=156;
numSub=22;
numSession=2;
numWinPerSubPerSeed=272;
numWinPerSeed=numWinPerSubPerSeed*numSub;
numSeed=4;
numWinPerSession=numWinPerSeed*numSeed;
numWin2Ses=numWinPerSession*numSession;
for i=1:numSession
    for j=1:numSeed
        startIndx=numWinPerSession*(i-1)+numWinPerSeed*(j-1)+numWinPerSubPerSeed*(subID-1)+1
        endIndx=numWinPerSession*(i-1)+numWinPerSeed*(j-1)+numWinPerSubPerSeed*subID
        b(startIndx:endIndx, :)=1000;
    end
end

bTransp=b';
bTransp(find(bTransp==1000))=[];
b=reshape(bTransp,numROI,[])';

zWinFullCorLasso21Sub=b;
size(zWinFullCorLasso21Sub)

save([resultDir, 'zWinFullCorLasso21Sub_', covType, '.mat'], 'zWinFullCorLasso21Sub')
