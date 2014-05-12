clear
clc
close all


session='session1'
if strcmp(session, 'session1')
    numState=5;
else
    numState=6;
end

resultDir=['/Users/zhenyang/Desktop/Zhen/results_4_3_13/all_10min/', session, filesep];

figDir='/Users/zhenyang/Desktop/Zhen/figs/';

%measure='transitions';
measure='TCPrctOverlapStateSum'
%measure='subSeedPrct';

% perform ttest between two seeds on the metric within a session
tmp=load([resultDir, measure, '_', session, '.mat']);
data=tmp.(measure);

pairs = nchoosek(1:size(data, 2), 2);
numPair=size(pairs, 1);

% Data in 3D matrix use the following code
pValue=zeros(numState,numPair)
for i=1:numPair
    for j=1:numState
        idx=pairs(i,:);
        x=squeeze(data(j, idx(1),:));
        y=squeeze(data(j, idx(2),:));
        [h,p] = ttest(x,y);
        pValue(j,i)=p;
    end
end
pValue

[pID,pN] = FDR(pValue,0.05)
negLogPValue=zeros(numState,numPair);
for i=1:numPair
    for j=1:numState
        if pValue(j,i) <= pID
            negLogPValue(j,i)=(-1)*log10(pValue(j,i));
        else
            negLogPValue(j,i)=0;
        end
    end
end

negLogPValue

% Data in 2D matrix use the following code
pValue=zeros(1,numPair)
for i=1:numPair
    idx=pairs(i,:);
    [h,p] = ttest(data(:, idx(1)), data(:, idx(2)));
    pValue(1,i)=p;
end
pValue

[pID,pN] = FDR(pValue,0.05)
negLogPValue=zeros(1,numPair);
for i=1:numPair
    if pValue(1,i) <= pID
        negLogPValue(1,i)=(-1)*log10(pValue(1,i));
    else
        negLogPValue(1,i)=0;
    end
end

negLogPValue
