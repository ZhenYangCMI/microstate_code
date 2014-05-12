clear
clc

subList={'3808535','8574662'};
numSub=length(subList);
numROI='0100';

for i=1:numSub
    sub=char(subList{i})
    mkdir (['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/testing/TR645_5min/', sub])
    resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/testing/TR645_5min/', sub, '/'];
    tmp=load(['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/testing/TR645_10min/', sub, '/ROISignals_k', numROI, '.mat']);
    ROISignal10min=tmp.ROISignals;
    ROISignals=ROISignal10min(1:442, :);
    size(ROISignals)
    save([resultDir, 'ROISignals_TR645k', numROI, '_5min.mat'], 'ROISignals')
end