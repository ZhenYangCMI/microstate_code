clear
clc
close all


numSub=21;
numCommonState=5;
numStateSes1=5;
numStateSes2=6;
numSession=2;

dataLength='all_10min';

resultDir='/Users/zhenyang/Desktop/Zhen/results/all_10min/2sesConcate/';
figDir='/Users/zhenyang/Desktop/Zhen/figs/5_24_13/';

%measureList={'transitions', 'TCPrctOverlapStateSum', 'prctDifNumStates'};
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
    tmp1=load([resultDir, measure,'_session1.mat']);
    data1=tmp1.(measure);
    meanData1=mean(data1);
    tmp2=load([resultDir, measure,'_session2.mat']);
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
    

    lineWidth=2
    figure(k)
    h=bar(ICCREML,0.6) % TCPrctOverlapStateSum used 0.6, other 2 measures used 0.4 to make them the same width
    set(gca,'LineWidth',lineWidth)
    %set(gca,'LineWidth',lineWidth,'xcolor','b','ycolor','b')
    set(h,'LineWidth',lineWidth,'Facecolor',[1,1,0],'EdgeColor',[128/255,128/255,128/255])
    if strcmp(measure, 'TCPrctOverlapStateSum')
        set(gca,'xTick',1:6, 'yTick', 0:0.2:1, 'box','off', 'xTickLabel',[],'yTickLabel',[]);
        xlim([0.5 6.5]);
        
    elseif strcmp(measure, 'prctDifNumStates')
        set(gca,'xTick',1:4, 'yTick', 0:0.2:1, 'box','off', 'xTickLabel',[],'yTickLabel',[]);
        xlim([0.5 4.5]);
        
        elseif strcmp(measure, 'transitions')
        set(gca,'xTick',1:4, 'yTick', 0:0.2:1, 'box','off', 'xTickLabel',[],'yTickLabel',[]);
        xlim([0.5 4.5]);
        
    end
    ylim([0 1]);
    %set(figure(1),'Colormap',mycmap), 
    %saveas(figure(k),[figDir,'ICC',measure,'2sesCmb21Sub.png'])
    
end

