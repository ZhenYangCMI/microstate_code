

clear
clc
close all

session='session1'
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;
numRand=100;

dataDir=['/home/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/', session, filesep];
resultsDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/',session,filesep];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/',session,filesep, 'surrogateData/']

numSeed=4;
numVol=884;
numWinPerSub=272;
numTotSub=22;
numROI=156;
numWinPerSeed=numWinPerSub*numTotSub;
TR='645';
tmp=load([dataDir,char(subID), '/ROISignals_seedROISignal.mat']);
allSeedTS=tmp.ROISignals;


randTSAllSeeds=zeros(numVol, numRand, numSeed);
for i=1:numSeed
    seedTS=squeeze(allSeedTS(:,i));
    [ out_ts_arr ] = autocorr_resample(seedTS, numRand );
    randTSAllSeeds(:,:,i)=out_ts_arr;
end

zwindFullCorLassoAllSurrogate=zeros(numWinPerSub*numSeed, numROI, numRand);
for j=1:numRand
    disp(['Working on rand ', num2str(j)])
    TC1=squeeze(randTSAllSeeds(:,j,:));
    lambdaList=load('/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/meanLambdaAtMaxLogAcrossWin.txt');
    lambdaSes1=squeeze(lambdaList(subOrder,1))
    disp('run the FCWinCreation Function')
    [winFullCor, winFullCorLasso, winPartialCorLasso]=FCWinCreationSurrogateDataSingleSub(TR,session,subID, lambdaSes1, TC1);
    zwindFullCorLassoAllSurrogate(:,:,j)=winFullCorLasso;
end
save([resultsDir, sprintf('zwindFullCorLassoAllSurrogate_sub%s_rand%d.mat', char(subID), numRand)], 'zwindFullCorLassoAllSurrogate')




