clear
clc
close all

session='session1'
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;
numRand=100;
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
dataDir=['/home/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/', session, filesep];
resultsDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/',session,filesep];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/',session,filesep, 'surrogateData/']

numVol=884;
numSeed=4;
numSub=22;
numROI=156;
TR='645';

corAllSeedAllSub=zeros(numSeed, numRand, numSub);
pAllSeedAllSub=zeros(numSeed, numRand, numSub);

for k=1:numSub
    sub=subList{k};
    disp (['Working on sub ', char(sub),' ......'])
    subDir=[analyDir,'data/',TR,'/all_10min/',session,'/',char(sub)];
    seedROISignals = load([subDir,'/ROISignals_seedROISignal.mat']);
    TC1tmp=seedROISignals.ROISignals;
        
    TC1=zeros(numVol, numRand, numSeed);
    for i=1:numSeed
        seedTS=squeeze(TC1tmp(:,i));
        [ out_ts_arr ] = autocorr_resample(seedTS, numRand );
        TC1(:,:,i)=out_ts_arr;
        for j=1:numRand
        [corr, pValue]=corrcoef(seedTS, out_ts_arr(:,j))
        corAllSeedAllSub(i, j, k)=corr(1,2);
        pAllSeedAllSub(i,j,k)=pValue(1,2);
        end
           end
    disp('correlations computed')
end
corAllSeedAllSub1D=reshape(corAllSeedAllSub, 1, []);
pAllSeedAllSub1D=reshape(pAllSeedAllSub, 1,[]);
[pID,pN] = FDR(pAllSeedAllSub1D,0.05)

numSigP=length(find(pAllSeedAllSub1D<pID));
threshold=corAllSeedAllSub(find(pAllSeedAllSub==pID))

figure(1)
hist(corAllSeedAllSub1D, -0.2:0.01:0.2)
title('Distribution of correlations between real TS and surrogated TS')
saveas(figure(1), [figDir, 'distribution of TScorrelation_realandsurrogate.png'])
