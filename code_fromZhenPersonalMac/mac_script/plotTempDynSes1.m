clear
clc
close all

session='session1';
numSub=22;
numState=5;
numSeed=4;
numSeedPair=6;

dataLength='all_10min';


resultDir='/Users/zhenyang/Desktop/Zhen/results/all_10min/session1/';
figDir='/Users/zhenyang/Desktop/Zhen/figs/5_24_13/';

% plot mean of 3D metrics
%measureList={'subSeedPrct', 'stateFreq', 'duration', 'TCPrctOverlapEachState'};
measureList={'duration'}
numMeasure=length(measureList);
colors=[0,100/255,0; 34/255,139/255,34/255;50/255,205/255,50/255;152/255,251/255,152/255;240/255,255/255,240/255;1,1,1];
colorMap = colorRamp(colors, 32)

for k=1:numMeasure
    measure=char(measureList{k});
    if strcmp(measure,'TCPrctOverlapEachState')
        numSeed=6;
    else
        numSeed=4;
    end
    % compute ICC and perform ttest on % time spent in a state
    tmp=load([resultDir, measure,'_session1.mat']);
    data=tmp.(measure);
    
    dataReshape=reshape(data, [], numSub)';
    meanReshap=mean(dataReshape);
    mean=reshape(meanReshap, numState, numSeed);
    
    % the state number needs to be reordered to match the 2 ses combined.
    % Here the order 1-5 is corresponding to the orginal: 4,2,5,1,3. The
    % new order 1-5 is corresponding to the original: 2,1,3,4,5 
    meanRecode=zeros(numState,numSeed);
        meanRecode(1,:)=mean(2,:);
        meanRecode(2,:)=mean(4,:);
        meanRecode(3,:)=mean(5,:);
        meanRecode(4,:)=mean(1,:);
        meanRecode(5,:)=mean(3,:);
    
    figure(k)
    imagesc(meanRecode)
    
    if strcmp(measure, 'stateFreq')
        caxis([0 4])
%         hcb=colorbar
%         set(hcb,'YTick',[0,2,4])
    elseif strcmp(measure, 'subSeedPrct')
        caxis([0 0.5])
%         hcb=colorbar
%         set(hcb,'YTick',[0,0.25,0.5])
    elseif strcmp(measure, 'duration') || strcmp(measure, 'TCPrctOverlapEachState')
        caxis([0 0.2])
%         hcb=colorbar
%         set(hcb,'YTick',[0,0.1,0.2])
    end
    if strcmp(measure, 'TCPrctOverlapEachState')
        %xlabel('Seed Pairs')
        set(gca,'xTick',0.5:6.5, 'yTick', 0.5:5.5,'xTickLabel',[],'yTickLabel',[]);
        xlim([0.5 6.5]);
    else
        %xlabel('Seeds')
        set(gca,'xTick',0.5:4.5, 'yTick', 0.5:4.5, 'xTickLabel',[],'yTickLabel',[]);
        xlim([0.5 4.5]);
    end
    ylim([0.5 5.5]);
    grid on
    set(gca,'GridLineStyle','-')
    set(gca,'LineWidth',2)
    set(figure(k),'Colormap',colorMap)
    %saveas(figure(k),[figDir,'Session1', measure,'.png'])
end




% plot boxplot for the 2D metrics

%measureList={'transitions', 'TCPrctOverlapStateSum', 'prctDifNumStates'};
measureList={'transitions'}
numMeasure=length(measureList);

for k=1:numMeasure
    measure=char(measureList{k});
    if strcmp(measure,'TCPrctOverlapStateSum')
        numSeed=6; %numSeed=numSeedPairs
    else
        numSeed=4;
    end
    % compute ICC and perform ttest on % time spent in a state
    tmp=load([resultDir, measure,'_',session,'.mat']);
    data=tmp.(measure);
    meanData1=mean(data);
    lineWidth=2
    figure(k)
    aboxplot(data,'labels',[5,10,15])
    %http://alex.bikfalvi.com/research/advanced_matlab_boxplot/
    % boxWidth: TCPrctOverlapStateSum 0.8, others 0.6
    %h=boxplot(data, 'boxstyle','outline','color',[0/255,128/255,0], 'medianstyle','line', 'widths',0.5)
    set(gca,'LineWidth',lineWidth)
    if strcmp(measure,'TCPrctOverlapStateSum')
        ylim([-0.2 1]);
        set(gca, 'yTick', -0.2:0.4:1,'xTickLabel',[],'yTickLabel',[], 'box','off');
        set(gca,'LineWidth',lineWidth)
        %set(h,'LineWidth',lineWidth)
        %     title(['session1'])
        %     ylabel('Pairs of seeds: TCPrctOverlapStateSum')
        %     xlabel('Seed Pairs')
    elseif strcmp (measure, 'transitions')
        ylim([-3 18]);
        set(gca, 'yTick', -3:6:18, 'xTickLabel',[],'yTickLabel',[], 'box','off');
        set(gca,'LineWidth',lineWidth)
        %set(h,'LineWidth',lineWidth)
        %     ylabel('Total number of transitions')
        %     xlabel('Seed')
    else
        ylim([-0.2 1]);
        set(gca,'yTick', -0.2:0.4:1, 'xTickLabel',[],'yTickLabel',[], 'box','off');
        set(gca,'LineWidth',lineWidth)
        %set(h,'LineWidth',lineWidth)
        %title(['session1'])
        %     ylabel('Prct time 4 seeds at different number of states ')
        %     xlabel('Number of states 4 seeds in at the same time')
    end
    %saveas(figure(k),[figDir, 'boxplot_',session,'_',measure,'.png'])
end


