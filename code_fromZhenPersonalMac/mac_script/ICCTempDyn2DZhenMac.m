clear
clc
close all


numSub=22;
numCommonState=5;
numStateSes1=5;
numStateSes2=6;
numSession=2;

dataLength='all_10min';
load('MyColormaps2','mycmap2') % for caxis >1
load('MyColormaps','mycmap') % for caxis [0 1]


resultDir1='/Users/zhenyang/Desktop/Zhen/results/all_10min/session1/';
resultDir2='/Users/zhenyang/Desktop/Zhen/results/all_10min/session2/';
figDir='/Users/zhenyang/Desktop/Zhen/figs/';

%measureList={'transitions', 'TCPrctOverlapStateSum'};
measureList={'prctDifNumStates'}
numMeasure=length(measureList);

for k=1:numMeasure
    measure=char(measureList{k});
    if strcmp(measure,'TCPrctOverlapStateSum')
        numSeed=6; %numSeed=numSeedPairs
    else
        numSeed=4;
    end
    % compute ICC and perform ttest on % time spent in a state
    tmp1=load([resultDir1, measure,'_session1.mat']);
    data1=tmp1.(measure);
    meanData1=mean(data1);
    tmp2=load([resultDir2, measure,'_session2.mat']);
    data2=tmp2.(measure);
    meanData2=mean(data2);
    
    ICCREML=zeros(1, numSeed);
    ICCIPN=zeros(1, numSeed);
    pValue=zeros(1,numSeed);
    for i=1:numSeed
        disp(['Work on seed ', num2str(i)])
        session1Data=squeeze(data1(:,i));
        session2Data=squeeze(data2(:,i));
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
    pValue
    
    [pID,pN] = FDR(pValue,0.05)
    negLogPValue=zeros(1, numSeed);
    for i=1:numSeed
        j=1;
        if pValue(j,i) <= pID
            if mean1(j,i)-mean2(j,i)>0;
                negLogPValue(j,i)=(-1)*log10(pValue(j,i));
            else
                negLogPValue(j,i)=log10(pValue(j,i));
            end
        else
            negLogPValue(j,i)=0;
        end
        
    end
    negLogPValue
end

% scatterplot all seeds/seedPairs and all states
if strcmp(measure,'TCPrctOverlapStateSum')
    numSeed=6;
else
    numSeed=4;
end

numICC=length(find(ICCREML>-0.01))
indx=find(ICCREML>-0.01);

for j=1:numICC
    figure(1)
    subplot(2,2,j)
    scatter(data1(:,j), data2(:,j))
    if strcmp(measure,'TCPrctOverlapStateSum')
        title(['SeedPair', num2str(indx(j))'])
        xlim([-0.5 1.5]);
        ylim([-0.5 1.5]);
        set(gca,'xTick',-0.5:0.5:1.5, 'yTick', -0.5:0.5:1.5);
    elseif strcmp (measure, 'transitions')
        title(['Seed', num2str(indx(j))'])
        xlim([-3 18]);
        ylim([-3 18]);
        set(gca,'xTick',-3:6:18, 'yTick', -3:6:18);
    else
        title(['4 seeds in ', num2str(indx(j)), 'state(s)'])
        xlim([-0.2 1]);
        ylim([-0.2 1]);
        set(gca,'xTick',-0.2:0.4:1, 'yTick', -0.2:0.4:1);
    end
    ylabel('Session 2')
    xlabel('Session 1')
        lsline
end
saveas(figure(1),[figDir,'scatter2SesAllSeedAllState_', measure, '.png'])



% scatterplot the data with ICCREML > 0.4
if strcmp(measure,'TCPrctOverlapStateSum')
    numICC=length(find(ICCREML>0.4))
    indx=find(ICCREML>0.4);
    if numICC~=0
        figure(1)
        for j=1:numICC
            subplot(1,3,j)
            scatter(data1(:,indx(j)), data2(:,indx(j)))
            xlim([-0.5 1.5]);
            ylim([-0.5 1.5]);
            set(gca,'xTick',-0.5:0.5:1.5, 'yTick', -0.5:0.5:1.5);
            ylabel('Session 2')
            xlabel('Session 1')
            title(['SeedPair', num2str(indx(j))'])
            lsline
        end
        saveas(figure(1),[figDir,'scatter2Sessions_', measure,'.png'])
    end
end


for k=1:numMeasure
    measure=char(measureList{k});
    figure(1)
    imagesc(data1)
    colorbar
    title(['Session1:', measure])
    ylabel('Subjects')
    
    if strcmp(measure, 'TCPrctOverlapStateSum')
        xlabel('Seed Pairs')
        set(gca,'xTick',1:6, 'yTick', 1:22);
        xlim([0.5 6.5]);
        caxis([0 1])
    else
        xlabel('Seeds')
        set(gca,'xTick',1:4, 'yTick', 1:22);
        xlim([0.5 4.5]);
        caxis([0 15])
    end
    set(figure(1),'Colormap',mycmap2)
    saveas(figure(1),[figDir,'Session1',measure,'.png'])
    
    figure(2)
    imagesc(data2)
    colorbar
    title(['Session2:', measure])
    ylabel('Subjects')
    caxis([0 15])
    if strcmp(measure, 'TCPrctOverlapStateSum')
        xlabel('Seed Pairs')
        set(gca,'xTick',1:6, 'yTick', 1:22);
        xlim([0.5 6.5]);
        caxis([0 1])
    else
        xlabel('Seeds')
        set(gca,'xTick',1:4, 'yTick', 1:22);
        xlim([0.5 4.5]);
        caxis([0 15])
    end
    set(figure(2),'Colormap',mycmap2)
    saveas(figure(2),[figDir,'Session2',measure,'.png'])
    
    figure(3)
    imagesc(ICCREML)
    colorbar
    title(['ICC',measure])
    caxis([0 1])
    if strcmp(measure, 'TCPrctOverlapStateSum')
        xlabel('Seed Pairs')
        set(gca,'xTick',1:6, 'yTick', []);
        xlim([0.5 6.5]);
    else
        xlabel('Seeds')
        set(gca,'xTick',1:4, 'yTick', []);
        xlim([0.5 4.5]);
    end
    ylim([0.5 1.5]);
    set(figure(3),'Colormap',mycmap)
    saveas(figure(3),[figDir,'ICC',measure,'.png'])
    
    figure(4)
    imagesc(negLogPValue)
    colorbar
    title(['Neglog10 of FDR corrected P Value: session dif on ', measure])
    caxis([-5 5])
    if strcmp(measure, 'TCPrctOverlapStateSum')
        xlabel('Seed Pairs')
        set(gca,'xTick',1:6, 'yTick', []);
        xlim([0.5 6.5]);
    else
        xlabel('Seeds')
        set(gca,'xTick',1:4, 'yTick', []);
        xlim([0.5 4.5]);
    end
    ylim([0.5 1.5]);
    saveas(figure(4),[figDir,'NegLogPValue', measure,'.png'])
    
    for m=1:numSession
        figNum=4+m;
        figure(figNum)
        data2Plot=eval(['meanData', num2str(m)]);
        imagesc(data2Plot)
        colorbar
        title(['Session',num2str(m), ': mean',measure])
        
        if strcmp(measure, 'TCPrctOverlapStateSum')
            xlabel('Seed Pairs')
            set(gca,'xTick',1:6, 'yTick', []);
            xlim([0.5 6.5]);
            caxis([0 1])
            set(figure(figNum),'Colormap',mycmap)
        else
            xlabel('Seeds')
            set(gca,'xTick',1:4, 'yTick', []);
            xlim([0.5 4.5]);
            caxis([0 15])
            set(figure(figNum),'Colormap',mycmap2)
        end
        ylim([0.5 1.5]);
        
        saveas(figure(figNum),[figDir,'session',num2str(m), '_mean',measure,'.png'])
    end
end
