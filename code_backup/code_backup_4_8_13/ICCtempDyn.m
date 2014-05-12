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

% compute ICC and perform ttest on % time spent in a state
tmp1=load([resultDir1, 'subSeedPrct_session1.mat']);
data1=tmp1.subSeedPrct;
tmp2=load([resultDir2, 'subSeedPrct_session2.mat']);
data2=tmp2.subSeedPrct;

ICCREMLSubSeedPrct=zeros(numCommonState,numSeed);
ICCIPNSubSeedPrct=zeros(numCommonState,numSeed);
for i=1:numSeed
    disp(['Work on seed ', num2str(i)])
    for j=1:numCommonState
        session1Data=squeeze(data1(j,i,:));
        session2Data=squeeze(data2(j,i,:));
        data=[session1Data,session2Data];
        Y=[session1Data;session2Data];
        time = [ones(numSub,1);2*ones(numSub,1)];
        sID=[[1:numSub]';[1:numSub]'];
        [ ICC1, idx_fail] = do_ICC(Y, time, [], [], sID);
        ICCREMLSubSeedPrct(j,i)=ICC1;
        ICC2 = IPN_icc(data,1,'single');
        ICCIPNSubSeedPrct(j,i)=ICC2;
        [h,p]=ttest(session1Data, session2Data);
        pValue(j,i)=p;
    end
end
ICCREMLSubSeedPrct
ICCIPNSubSeedPrct
data1Reshape=reshape(data1, [], numSub)';
mean1=mean(data1Reshape);
mean1SubSeedPrct=reshape(mean1, numStateSes1, numSeed);
data2Reshape=reshape(data2, [], numSub)';
mean2=mean(data2Reshape);
mean2SubSeedPrct=reshape(mean2, numStateSes2, numSeed);

[pID,pN] = FDR(pValue,0.05)
for i=1:numSeed
    for j=1:numCommonState
                if pValue(j,i) <= pID
                    negLogPValue(j,i)=(-1)*log10(pValue(j,i));
                    i           if 
        else
            pValue(j,i)=0;
        end
    end
end


figure(1)
imagesc(ICCREMLSubSeedPrct)
colorbar
title('ICC %time spent in a state')
xlabel('Seeds')
ylabel('States')
caxis([0 1])
set(gca,'xTick',1:4, 'yTick', 1:5);
xlim([0.5 4.5]);
ylim([0.5 5.5]);
set(figure(1),'Colormap',mycmap)
saveas(figure(1),[figDir,'ICCSubSeedPrct.png'])

figure(2)
imagesc(mean1SubSeedPrct)
colorbar
title('Session 1: %time spent in a state')
xlabel('Seeds')
ylabel('States')
caxis([0 1])
set(gca,'xTick',1:4, 'yTick', 1:5);
xlim([0.5 4.5]);
ylim([0.5 5.5]);
set(figure(2),'Colormap',mycmap)
saveas(figure(2),[figDir,'Session1SubSeedPrct.png'])

figure(3)
imagesc(mean2SubSeedPrct)
colorbar
title('Session 2: %time spent in a state')
xlabel('Seeds')
ylabel('States')
caxis([0 1])
set(gca,'xTick',1:4, 'yTick', 1:6);
xlim([0.5 4.5]);
ylim([0.5 6.5]);
set(figure(3),'Colormap',mycmap)
saveas(figure(3),[figDir,'Session2SubSeedPrct6States.png'])


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




