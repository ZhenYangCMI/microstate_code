clear
clc
close all

covType='compCor'
session='session1'
dataLength='all_10min';
winSize = 69;
% window sliding step in TRs
step=3;
% total number of time points
N_vol=884;
numSub=22;
numSeed=4;
% number of windows
numWinPerSub=floor((N_vol-winSize)/step)+1;
numWinPerSeed=numWinPerSub*numSub
Time=ceil(winSize/2):step:(ceil(winSize/2)+step*(numWinPerSub-1));
% seedPairs=seed1 and seed2; 1&3; 1&4; 2&3; 2&4; 3&4 
numState=4;

numROI=156;
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,'/', covType, '/', session,'/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, '/', covType, '/', session, '/'];

indx=load([resultDir,'clustIndxNormWinAllSeeds_FullCorLasso_',session,'win', num2str(winSize),'_', covType, '.txt']);
numClust=length(unique(indx));
disp ('Files loaded successfully.')


% compute the %time spent in each state for each seed and each sub
subSeedPrct=zeros(numState, numSeed, numSub);
for n = 1:numState
    for i = 1:numSeed
        for j = 1:numSub
            begin_ndx=(i-1)*numSub*numWinPerSub+(j-1)*numWinPerSub+1;
            end_ndx=begin_ndx+numWinPerSub-1;
            subSeed(n,i,j)=sum(indx(begin_ndx:end_ndx)==n);
            subSeedPrct(n,i,j)=sum(indx(begin_ndx:end_ndx)==n)/numWinPerSub*100;
            %disp(sprintf('begin %d, end %d, val %f',begin_ndx,end_ndx,subSeedPrct(n,i,j)));
        end
    end
end
save([resultDir, 'subSeedPrct_', session,'_', covType, '.mat'], 'subSeedPrct')

subSeedPrct2D=reshape(subSeedPrct, [], numSub)';
average=mean(subSeedPrct2D);
meanTotalDuration=reshape(average, numState, numSeed)

figure(1)
imagesc(meanTotalDuration)
colorbar
caxis([0 0.7])
title([covType])
set(gca,'xTick',1:4, 'yTick', 1:4);
xlabel('Seeds')
ylabel('states')

