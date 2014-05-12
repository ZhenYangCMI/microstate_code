
clear all
clc


% for windSize 136 only 1 sub was removed. For winsize34, we need to remove
% 2 subjects, the following code was run to rm one subject
winSize=136;
a=load(['zWinFullCorLasso_OptimalLambdaPerSub_645_2sessions_allsub_win', num2str(winSize), '.mat'])
b=a.zWinFullCorLasso;

% x=randsample(22,1)
% Get an answer of 18, subject 18 will be excluded from this analysis due
% to R object limit
% subID=x; The sub being removed
subID=18;
numROI=156;
numSub=22;
numSession=2;
if winSize==34
numWinPerSubPerSeed=284; 
elseif winSize==136
    numWinPerSubPerSeed=250;
end
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

save(['zWinFullCorLasso_OptimalLambdaPerSub_645_2sessions_21sub_win', num2str(winSize), '.mat'], 'zWinFullCorLasso21Sub')


% the following code was run to remove the second sub (subNum12: 3313349: who dosen't have NP data and had the highest motion among the subject who had NP)

clear
clc

winSize=34;
a=load(['zWinFullCorLasso_OptimalLambdaPerSub_645_2sessions_21sub_win', num2str(winSize), '.mat'])
b=a.zWinFullCorLasso21Sub;

% the sub being removed
subID=12;
numROI=156;
numSub=21;
numSession=2;
if winSize==34
numWinPerSubPerSeed=284; 
elseif winSize==136
    numWinPerSubPerSeed=250;
end
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

zWinFullCorLasso20Sub=b;
size(zWinFullCorLasso20Sub)

save(['zWinFullCorLasso_OptimalLambdaPerSub_645_2sessions_20sub_win', num2str(winSize), '.mat'], 'zWinFullCorLasso20Sub')