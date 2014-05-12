clear
clc
close all

numSeed=4;
numSub=22;
numCommonState=5;
numStateSes1=5;
numStateSes2=6;

dataLength='all_10min';
load('MyColormaps2','mycmap2')
load('MyColormaps','mycmap')


resultDir1=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,filesep, 'session1/'];
resultDir2=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,filesep, 'session2/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep, 'correlTwoSessions/'];

% compute ICC and perform t-test on totoal number of transition
tmp1=load([resultDir1, 'transitions_session1.mat']);
transit1=tmp1.transitions;
tmp2=load([resultDir2, 'transitions_session2.mat']);
transit2=tmp2.transitions;

figure(1)
imagesc(transit1)
colorbar
title('Number of transitions ')
xlabel('Seeds')
ylabel('Subjects')
caxis([0 15])
set(gca,'xTick',1:4);
xlim([0.5 4.5]);
set(figure(1),'Colormap',mycmap2)
saveas(figure(1),[figDir,'TotNumTransitSession1.png'])

figure(2)
imagesc(transit2)
colorbar
title('Number of transitions ')
xlabel('Seeds')
ylabel('Subjects')
caxis([0 15])
set(gca,'xTick',1:4);
xlim([0.5 4.5]);
set(figure(2),'Colormap',mycmap2)
saveas(figure(2),[figDir,'TotNumTransitSession2.png'])

ICCREML=zeros(1, numSeed);
ICCIPN=zeros(1, numSeed);
pValue=zeros(1,numSeed);
for i=1:numSeed
    disp(['Work on seed ', num2str(i)])
    session1Data=squeeze(transit1(:,i));
    session2Data=squeeze(transit2(:,i));
    data=[session1Data,session2Data];
    Y=[session1Data;session2Data];
    time = [ones(numSub,1);2*ones(numSub,1)];
    sID=[[1:numSub]';[1:numSub]'];
    [ ICC1, idx_fail] = do_ICC(Y, time, [], [], sID);
    ICCREML(1,i)=ICC1;
    ICC2 = IPN_icc(data,1,'single')
    ICCIPN(1,i)=ICC2;
    [h,p]=ttest(session1Data, session2Data);
    pValue(1,i)=p;
end
ICCREML
ICCIPN