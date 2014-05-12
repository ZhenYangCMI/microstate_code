clear
clc
close all


numSub=22;
numCommonState=5;
numStateSes1=5;
numStateSes2=6;
session='session1'

dataLength='all_10min';
load('MyColormaps2','mycmap2')
load('MyColormaps','mycmap')

resultDir=['/Users/zhenyang/Desktop/Zhen/results_4_3_13/all_10min/', session, filesep];

figDir='/Users/zhenyang/Desktop/Zhen/figs/';

measure1='transitions1';
measure2='transitions2';
measure3='transitions3';

% measure1='TCPrctOverlap3Period1';
% measure2='TCPrctOverlap3Period2';
% measure3='TCPrctOverlap3Period3';

% perform ANOVA on the metrics
tmp1=load([resultDir, measure1, '_', session, '.mat']);
data1=tmp1.(measure1);
tmp2=load([resultDir, measure2, '_', session, '.mat']);
data2=tmp2.(measure2);
tmp3=load([resultDir, measure3, '_', session, '.mat']);
data3=tmp3.(measure3);

if strcmp(measure1,'TCPrctOverlap3Period1')
    numSeed=6;
else
    numSeed=4;
end

% compute the ICC between 2 period within one session: period 1 and period
% 2, period 2 and period 3
ICCREML=zeros(1, numSeed);
for i=1:numSeed
    disp(['Work on seed ', num2str(i)])
    session1Data=squeeze(data1(:,i));
    session2Data=squeeze(data2(:,i));
    session3Data=squeeze(data3(:,i));
    Y1=[session1Data;session2Data];
    Y2=[session2Data;session3Data];
    time = [ones(numSub,1);2*ones(numSub,1)];
    sID=[[1:numSub]';[1:numSub]'];
    [ ICC1, idx_fail1] = do_ICC(Y1, time, [], [], sID);
    [ ICC2, idx_fail2] = do_ICC(Y2, time, [], [], sID);
    ICCREML1(1,i)=ICC1;
    ICCREML2(1,i)=ICC2;
end
ICCREML1
ICCREML2


numCommonState=1;
pValue=zeros(numCommonState, numSeed);
for i=1:numSeed
    disp(['Work on seed ', num2str(i)])
    for j=1:numCommonState
        period1=squeeze(data1(:,i));
        period2=squeeze(data2(:,i));
        period3=squeeze(data3(:,i));
        cat=[period1,period2,period3];
        [p,table,stats]=anova1(cat);
        pValue(j,i)=p;
    end
end
pValue

[pID,pN] = FDR(pValue,0.05)
negLogPValue=zeros(numCommonState, numSeed);
for i=1:numSeed
    for j=1:numCommonState
        if pValue(j,i) <= pID
            negLogPValue(j,i)=(-1)*log10(pValue(j,i));
        else
            negLogPValue(j,i)=0;
        end
    end
end
negLogPValue



% plot figures
if strcmp(session, 'session1')
    numState=5;
else
    numState=6;
end

mean1=mean(data1);
mean2=mean(data2);
mean3=mean(data3);

close all
numPeriod=3;

for i=1:numPeriod
    eval(['meanData=mean' num2str(i) ';']);
    eval(['measure=measure' num2str(i) ';'])
    figure(i)
    imagesc(meanData)
    colorbar
    title([measure])
    if strcmp(measure,'TCPrctOverlap3Period1') || strcmp(measure,'TCPrctOverlap3Period2')||strcmp(measure,'TCPrctOverlap3Period3')
        caxis([0 1])
        xlim([0.5 6.5]);
        set(gca,'xTick',1:6, 'yTick', []);
        xlabel('Seed Pairs')
    else
        caxis([0 4])
        xlim([0.5 4.5]);
        set(gca,'xTick',1:4, 'yTick', []);
        xlabel('Seeds')
    end
    set(figure(i),'Colormap',mycmap)
    saveas(figure(i),[figDir,session, 'mean', measure,'.png'])
end

figure(4)
imagesc(negLogPValue)
colorbar
title(['Neglog10 of FDR corrected P Value: periods dif on ', measure1(1:end-1)])
if strcmp(measure1,'TCPrctOverlap3Period1')
    caxis([0 1])
    xlim([0.5 6.5]);
    set(gca,'xTick',1:6, 'yTick', []);
    xlabel('Seed Pairs')
else
    caxis([0 4])
    xlim([0.5 4.5]);
    set(gca,'xTick',1:4, 'yTick', []);
    xlabel('Seeds')
end
figName=[session, '_NegLogPValue_', measure1(1:end-1), '.png'];
saveas(figure(4),[figDir, figName])


for i=1:numPeriod
    eval(['plotData=data' num2str(i) ';']);
    eval(['measure=measure' num2str(i) ';'])
    figure(4+i)
    
    imagesc(plotData)
    colorbar
    title([measure])
    if strcmp(measure,'TCPrctOverlap3Period1') || strcmp(measure,'TCPrctOverlap3Period2')||strcmp(measure,'TCPrctOverlap3Period3')
        caxis([0 1])
        xlim([0.5 6.5]);
        set(gca,'xTick',1:6, 'yTick', 1:2:21);
        xlabel('Seed Pairs')
    else
        caxis([0 4])
        xlim([0.5 4.5]);
        set(gca,'xTick',1:4, 'yTick', 1:2:21);
        xlabel('Seeds')
    end
    ylabel('Subjects')
    ylim([0.5 22.5])
    figNum=4+i;
    set(figure(figNum),'Colormap',mycmap)
    saveas(figure(figNum),[figDir,session, measure,'.png'])
end

close all
periodPair=[1,2;2,3]
numICC=2;
for i=1:numICC
    p1=squeeze(periodPair(i,1));
    p2=squeeze(periodPair(i,2))
    figure(8+i)
    plotData=eval(['ICCREML' num2str(i)])
    imagesc(plotData)
    colorbar
    title(['ICC between period ', num2str(p1), ' and period', num2str(p2)])
    if strcmp(measure1,'TCPrctOverlap3Period1')
        xlim([0.5 6.5]);
        set(gca,'xTick',1:6, 'yTick', []);
        xlabel('Seed Pairs')
    else
        xlim([0.5 4.5]);
        set(gca,'xTick',1:4, 'yTick', []);
        xlabel('Seeds')
    end
    caxis([0 1])
    figName=sprintf('%s_ICC%d%d_%s.png', session,p1, p2, measure1(1:end-1));
    figNum=8+i;
    set(figure(figNum),'Colormap',mycmap)
    saveas(figure(figNum),[figDir, figName])
end

% scatterplot the data with ICCREML > 0.4
for i=1:numICC
    ICC1=squeeze(periodPair(i,1));
    ICC2=squeeze(periodPair(i,2));
    ICCREML=eval(['ICCREML' num2str(i) ';']);
    numSigICC=length(find(ICCREML>0.4))
    indx=find(ICCREML>0.4);
    dataPlot1=eval(['data' num2str(ICC1)]);
    dataPlot2=eval(['data' num2str(ICC2)]);
    
    if numSigICC~=0
        figure(10+i)
        for j=1:numSigICC
            subplot(3,3,j)
            scatter(dataPlot1(:,indx(j)), dataPlot2(:,indx(j)))
                                    ylabel(['Period',num2str(ICC2)])
            xlabel(['Period',num2str(ICC1)])
            if strcmp(measure1,'TCPrctOverlap3Period1')
                title(['SeedPair', num2str(indx(j))'])
                xlim([-0.5 1.5]);
            ylim([-0.5 1.5]);
            set(gca,'xTick',-0.5:0.5:1.5, 'yTick', -0.5:0.5:1.5);
            else
                title(['Seed', num2str(indx(j))'])
                xlim([-2 8]);
            ylim([-2 8]);
            set(gca,'xTick',-2:2:8, 'yTick', -2:2:8);
            end
            lsline
        end
    end
    figNum=10+i;
    figName=sprintf('%s_scatterPeriod%d%d_%s.png', session,ICC1, ICC2, measure1(1:end-1));
    saveas(figure(figNum),[figDir,figName])
end


