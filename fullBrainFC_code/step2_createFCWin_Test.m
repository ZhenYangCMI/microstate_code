%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script do the following for each session separately:
%1. generate three types of Fisher Z transformed FC windows: Pearson's r, gLasso full, and partial covariance windows
%Each type of the FC windows from all subjects are concatenated into one file.
%2. extract the features from the full brain correlation matrix to form the
%featureWin and Z standardize each window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

% load subList
subList={'3808535','8574662'};
numSub=length(subList);
numROI='0200';
TR='645';
session='session1'

% define the analysis strategy, use to fid the right TS data
covType='GSR';

% define the windodw parameters
% below for TR645
winSize=69;
step=3; % in TR
numVol=442;

% winSize=18;
% step=1; % in TR
% numVol=116;

ROISignalName=['ROISignals_TR', TR, 'k', numROI, '_5min.mat'];

%define data and resultDir
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/'];
dataDir=[analyDir, 'testing/TR', TR, '_5min/'];
resultDir=[analyDir,'testing/results/TR', TR, '_5min/'];

% load the lambdaList
tmp=load([analyDir,'testing/results/TR', TR, '_5min/meanLambdaAtMaxLog_', numROI, '_', session, '_2subjects_TR', TR, '.mat']);
lambdaList=tmp.meanLambdaAtMaxLogAcrossWin

% 1 generate full and partial correlation windows
disp(['Create full and partial correlation windows for', session,'.'])
%run the FCWinCreation function to get the concatenated Fisher Z
% transformed FC windows
clear winFullCor winFullCorLasso winPartialCorLasso
[winFullCor, winFullCorLasso, winPartialCorLasso]=FCWinCreation(dataDir,subList, lambdaList, winSize, step, numVol, ROISignalName);

% save these FC windows
save([resultDir,'winFullCor_Pearson_',numROI, '_', session,'.mat'],'winFullCor', '-v7.3')
 save([resultDir,'winFullCorLasso_OptimalLambdaPerSub_',session,'win', num2str(winSize), '_',covType,'.mat'],'winFullCorLasso', '-v7.3')
save([resultDir,'winPartialCorLasso_OptimalLambdaPerSub_', session,'.mat'],'winPartialCorLasso', '-v7.3')


% 2 Extract the features from the full brain correlation matrix and
% standardize the features for each featureWin

%tmp1=load([resultDir,'winFullCor_Pearson_',numROI, '_', session,'.mat'])
%tmp2=load([resultDir,'winFullCorLasso_OptimalLambdaPerSub_',session,'win', num2str(winSize), '_',covType,'.mat'])
%tmp3=load([resultDir,'winPartialCorLasso_OptimalLambdaPerSub_', session,'.mat'])
%winFullCor=tmp1.winFullCor;
%winFullCorLasso=tmp2.winFullCorLasso;
%winPartialCorLasso=tmp3.winPartialCorLasso; 


winType={'winFullCor','winFullCorLasso', 'winPartialCorLasso'}; %'winFullCor', 'winFullCorLasso', 'winPartialCorLasso'

for i=1:length(winType)
    disp(['Working on ', char(winType{i})])
    FCWinFullBrain=eval(char(winType{i}));
    [featureWin zfeatureWin] = featureExtract( FCWinFullBrain, subList, winSize, step, numVol );
   save([resultDir,'zfeatureFC_', char(winType{i}), '_',numROI, '_', session,'.mat'],'zfeatureWin')
   save([resultDir,'featureFC_', char(winType{i}), '_', numROI, '_', session,'.mat'],'featureWin')
end
