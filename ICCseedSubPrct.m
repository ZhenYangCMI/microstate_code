clear
clc
close all

numSeed=4;
numSub=22;
% the first column in each row of statePair is the state number of session1
% the second colume in each row is the state number of session2
% each row is a pair of highly correlated states across two sessions
statePairList=[4,1;5,5;3,6;1,2;1,3;2,2;2,3];
numStatePairs=size(statePairList,1);

dataLength='all_10min';
load('MyColormaps','mycmap')


resultDir1=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,filesep, 'session1/'];
resultDir2=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,filesep, 'session2/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep, 'correlTwoSessions/'];

% session 1: compute percent of each state
indx1=load([resultDir1,'clustIndxNormWinAllSeeds_FullCorLasso_session1_10min.txt']);
numClust1=length(unique(indx1));
[seedPrct1, subPrct1, subSeedPrct1] = DoQuickStats(indx1,numClust1);

% session 2: compute percent of each state
indx2=load([resultDir2,'clustIndxNormWinAllSeeds_FullCorLasso_session2_10min.txt']);
numClust2=length(unique(indx2));
[seedPrct2, subPrct2, subSeedPrct2] = DoQuickStats(indx2,numClust2);

ICCoutput=zeros(numSeed,numStatePairs);
for i=1:numSeed
    disp(['Work on seed ', num2str(i)])
    for j=1:numStatePairs
        session1StateNum=statePairList(j,1);
        session2StateNum=statePairList(j,2);
        session1Data=squeeze(subSeedPrct1(session1StateNum,i,:));
        session2Data=squeeze(subSeedPrct2(session2StateNum,i,:));
        data=[session1Data,session2Data];
        Y=[session1Data;session2Data];
        time = [ones(numSub,1);2*ones(numSub,1)];
        sID=[[1:numSub]';[1:numSub]'];
        [ ICC1, idx_fail] = do_ICC(Y, time, [], [], sID);
        ICCoutput1(i,j)=ICC1
        ICC2 = IPN_icc(data,3,'k')
        ICCoutput2(i,j)=ICC2
    end
end

figure(1)
imagesc(ICCoutput1)
colorbar
title('ICC between sessions on the percent windows in a state (ICCoutput1)')
xlabel('States exhibiting correlation between sessions')
ylabel('Seeds')
caxis([-1 1])
set(gca,'YTick',1:4);
ylim([0.5 4.5]);
saveas(figure(1),[figDir,'ICC_REML.png'])

figure(2)
imagesc(ICCoutput2)
colorbar
title('ICC between sessions on the percent windows in a state (ICCoutput2)')
xlabel('States exhibiting correlation between sessions')
ylabel('Seeds')
caxis([-1 1])
set(gca,'YTick',1:4);
ylim([0.5 4.5]);
saveas(figure(2),[figDir,'ICC_IPN.png'])

