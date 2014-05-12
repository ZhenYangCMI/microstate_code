clear
clc
close all


numSub=22;
numCommonState=5;
numStateSes1=5;
numStateSes2=6;
session='session1'
numSeed=4;

dataLength='all_10min';
load('MyColormaps2','mycmap2')
load('MyColormaps','mycmap')

resultDir=['/Users/zhenyang/Desktop/Zhen/results/all_10min/'];
resultDir1=['/Users/zhenyang/Desktop/Zhen/results/all_10min/session1/'];
resultDir2=['/Users/zhenyang/Desktop/Zhen/results/all_10min/session2/'];

figDir='/Users/zhenyang/Desktop/Zhen/figs/4_28_13';


% perform ANOVA on the metrics
tmp1=load([resultDir1, 'transitions_session1.mat']);
transitions1=tmp1.transitions;
tmp2=load([resultDir2, 'transitions_session2.mat']);
transitions2=tmp2.transitions;
save([resultDir1, 'transitions.session1.txt'],'transitions1')
save([resultDir2, 'transitions.session2.txt'],'transitions2')
transitions=[transitions1, transitions2];
save([resultDir, 'transitions.2sessions.txt'],'-ascii', '-double', '-tabs', 'transitions')

pValue=zeros(numCommonState, numSeed);
for i=1:numSeed
    disp(['Work on seed ', num2str(i)])
    for j=1:numCommonState
        period1=squeeze(data1(j,i,:));
        period2=squeeze(data2(j,i,:));
        period3=squeeze(data3(j,i,:));
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

if strcmp(session, 'session1')
    numState=5;
else
    numState=6;
end
data1Reshape=reshape(data1, [], numSub)';
mean1Reshap=mean(data1Reshape);
mean1=reshape(mean1Reshap, numState, numSeed);
data2Reshape=reshape(data2, [], numSub)';
mean2Reshap=mean(data2Reshape);
mean2=reshape(mean2Reshap, numState, numSeed);
data3Reshape=reshape(data3, [], numSub)';
mean3Reshap=mean(data3Reshape);
mean3=reshape(mean3Reshap, numState, numSeed);

close all
numPeriod=3;

for i=1:numPeriod
    eval(['meanData=mean' num2str(i) ';']);
    eval(['measure=measure' num2str(i) ';'])
    figure(i)
    imagesc(meanData)
    colorbar
    title([measure])
    ylabel('States')
    caxis([0 1])
    xlim([0.5 4.5]);
    xlabel('Seeds')
    if strcmp(session, 'session1')
        set(gca,'xTick',1:4, 'yTick', 1:5);
        ylim([0.5 5.5]);
    else
        set(gca,'xTick',1:4, 'yTick', 1:6);
        ylim([0.5 6.5]);
    end
    set(figure(i),'Colormap',mycmap)
    saveas(figure(i),[figDir,session, measure,'.png'])
end

figure(4)
imagesc(negLogPValue)
colorbar
title(['Neglog10 of FDR corrected P Value: periods dif on subSeedPrct'])
ylabel('States')
caxis([-5 5])
xlim([0.5 4.5]);
xlabel('Seeds')
set(gca,'xTick',1:4, 'yTick', 1:5);
ylim([0.5 5.5]);
saveas(figure(4),[figDir,session, '_NegLogPValue on subSeedPrct.png'])


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
    p1=squeeze(periodPair(i,1))
    p2=squeeze(periodPair(i,2))
    ICCREML=eval(['ICCREML' num2str(i) ';']);
    numSigICC=length(find(ICCREML>0.6))
    indx=find(ICCREML>0.6);
    dataPlot1=eval(['data' num2str(p1)]);
    dataPlot2=eval(['data' num2str(p2)]);
    
    if numSigICC~=0
        figure(10+i)
        for j=1:numSigICC
            subplot(4,3,j)
            col=ceil(indx(j)/numCommonState)
            if rem(indx(j), numCommonState)~=0
                row=rem(indx(j), numCommonState)
            else
                row=numCommonState
            end
            scatter(dataPlot1(row,col,:), dataPlot2(row,col,:))
            ylabel(['Period',num2str(p2)])
            xlabel(['Period',num2str(p1)])
            title(['Seed', num2str(col), 'state', num2str(row)])
            xlim([-0.5 1.5]);
            ylim([-0.5 1.5]);
            set(gca,'xTick',-0.5:0.5:1.5, 'yTick', -0.5:0.5:1.5);
            lsline
        end
    end
    figNum=10+i;
    figName=sprintf('%s_scatterPeriod%d%d_%s.png', session,p1, p2, measure1(1:end-1));
    saveas(figure(figNum),[figDir,figName])
end






