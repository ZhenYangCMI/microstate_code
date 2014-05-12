

clear
clc
close all

session='session1'
subID={'8574662'} % the subject with lowest ('8574662') and highest ('3795193') motion was used to do the pilot analysis
subOrder=20; % subOrder for '8574662' is 20, for '3795193' is 14.
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
    [winFullCor, winFullCorLasso, winPartialCorLasso]=FCWinCreationSurrogateData(TR,session,subID, lambdaSes1, TC1);
    zwindFullCorLassoAllSurrogate(:,:,j)=winFullCorLasso;
end
save([resultsDir, sprintf('zwindFullCorLassoAllSurrogate_sub%s_rand%d.mat', char(subID), numRand)], 'zwindFullCorLassoAllSurrogate')

%tmp1=load([resultsDir, 'zwindFullCorLassoAllSurrogate_sub', char(subID), '_rand', num2str(numRand), '.mat']);
%zwindFullCorLassoAllSurrogate=tmp1.zwindFullCorLassoAllSurrogate;
FCwin=load([resultsDir, 'zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'.mat'])
FCwin=FCwin.zWinFullCorLasso;

%close all
for i=1:numSeed
    original1Seed=FCwin(1+numWinPerSub*(subOrder-1)+numWinPerSeed*(i-1):numWinPerSub*subOrder+numWinPerSeed*(i-1), :);
    original1Seed1D=reshape(original1Seed, [],1);
    surrogate1Seed=zwindFullCorLassoAllSurrogate(1+numWinPerSub*(i-1):numWinPerSub*i, :, :);
    surrogate1Seed1D=reshape(surrogate1Seed, [],1);
    
    figure(i)
        subplot(1,2,1)
     
    boxplot(original1Seed1D)
    ylim([-5, 5])
       set(gca, 'xTick',[])
    xlabel('Original Data')
    ylabel('Fisher z scores')
    title(['seed ', num2str(i)])
    subplot(1,2,2)
    boxplot(surrogate1Seed1D)
     ylim([-5, 5])
     xlabel(['Surrogated Data (N=', num2str(numRand), ')'])
     title(['seed ', num2str(i)])
    saveas(figure(i),[figDir, sprintf('var_orginalVSsurrogate_sub%s_seed%d_rand%d.png', char(subID), i, numRand)])
end

