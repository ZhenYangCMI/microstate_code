clear
clc
close all


numSub=22;
numCommonState=5;
numStateSes1=5;
numStateSes2=6;
numSession=2;
numSeed=4;
numSeedPair=6;

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
        numLoop=numSeedPair; %numSeed=numSeedPairs
    else
        numLoop=numSeed;
    end
    % compute ICC and perform ttest on % time spent in a state
    tmp1=load([resultDir1, measure,'_session1.mat']);
    data1=tmp1.(measure);
    meanData1=mean(data1);
    tmp2=load([resultDir2, measure,'_session2.mat']);
    data2=tmp2.(measure);
    meanData2=mean(data2);
    
    ICCREML=zeros(1, numLoop);
    ICCIPN=zeros(1, numLoop);
    pValue=zeros(1,numLoop);
    for i=1:numLoop
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
    negLogPValue=zeros(1, numLoop);
    for i=1:numLoop
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


for k=1:numMeasure
    measure=char(measureList{k});
    lineWidth=2
    figure(1)
    h=bar(ICCREML)
    set(h,'LineWidth',lineWidth)
    set(gca,'LineWidth',lineWidth, 'box','off')
    if strcmp(measure, 'TCPrctOverlapStateSum')
        xlabel('Seed Pairs')
        set(gca,'xTick',1:6, 'yTick', 0:0.2:1);
        xlim([0.5 6.5]);
    else
        xlabel('Seeds')
        set(gca,'xTick',1:4, 'yTick', 0:0.2:1);
        xlim([0.5 4.5]);
    end
    ylim([0 1]);
    set(figure(3),'Colormap',mycmap)
    saveas(figure(1),[figDir,'ICC',measure,'.png'])
    
    for m=1:numSession
        figNum=1+m;
        data2Plot=eval(['data', num2str(m)]);
                figure(figNum)
        h=boxplot(data2Plot)
        if strcmp(measure,'TCPrctOverlapStateSum')
            ylim([-0.2 1]);
            set(gca, 'yTick', -0.2:0.4:1,'xTickLabel',[],'yTickLabel',[], 'box','off');
            set(gca,'LineWidth',lineWidth)
            set(h,'LineWidth',lineWidth)
            %     title(['session1'])
            %     ylabel('Pairs of seeds: TCPrctOverlapStateSum')
            %     xlabel('Seed Pairs')
        elseif strcmp (measure, 'transitions')
            ylim([-3 18]);
            set(gca, 'yTick', -3:6:18, 'xTickLabel',[],'yTickLabel',[], 'box','off');
            set(gca,'LineWidth',lineWidth)
            set(h,'LineWidth',lineWidth)
            %     ylabel('Total number of transitions')
            %     xlabel('Seed')
        else
            ylim([-0.2 1]);
            set(gca,'yTick', -0.2:0.4:1, 'xTickLabel',[],'yTickLabel',[], 'box','off');
            set(gca,'LineWidth',lineWidth)
            set(h,'LineWidth',lineWidth)
            %title(['session1'])
            %     ylabel('Prct time 4 seeds at different number of states ')
            %     xlabel('Number of states 4 seeds in at the same time')
        end
        
        saveas(figure(figNum),[figDir,'session',num2str(m), '_mean',measure,'.png'])
    end
    
    
end
