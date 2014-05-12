clear
clc
close all


numSub=22;
numCommonState=5;
numStateSes1=5;
numStateSes2=6;

dataLength='all_10min';
load('MyColormaps2','mycmap2')
load('MyColormaps','mycmap')


resultDir1='/Users/zhenyang/Desktop/Zhen/results_4_3_13/all_10min/session1/';
resultDir2='/Users/zhenyang/Desktop/Zhen/results_4_3_13/all_10min/session2/';
figDir='/Users/zhenyang/Desktop/Zhen/figs/';

%measureList={'subSeedPrct', 'stateFreq', 'duration', 'TCPrctOverlapEachState'};
measureList={'duration'}
numMeasure=length(measureList);


for k=1:numMeasure
    measure=char(measureList{k});
    if strcmp(measure,'TCPrctOverlapEachState')
        numSeed=6;
    else
        numSeed=4;
    end
    % compute ICC and perform ttest on % time spent in a state
    tmp1=load([resultDir1, measure,'_session1.mat']);
    data1=tmp1.(measure);
    tmp2=load([resultDir2, measure,'_session2.mat']);
    data2=tmp2.(measure);
    
    ICCREML=zeros(numCommonState,numSeed);
    ICCIPN=zeros(numCommonState,numSeed);
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
            ICCREML(j,i)=ICC1;
            ICC2 = IPN_icc(data,1,'single');
            ICCIPN(j,i)=ICC2;
            [h,p]=ttest(session1Data, session2Data);
            pValue(j,i)=p;
        end
    end
    ICCREML
    ICCIPN
    pValue
    
    data1Reshape=reshape(data1, [], numSub)';
    mean1Reshap=mean(data1Reshape);
    mean1=reshape(mean1Reshap, numStateSes1, numSeed);
    data2Reshape=reshape(data2, [], numSub)';
    mean2Reshap=mean(data2Reshape);
    mean2=reshape(mean2Reshap, numStateSes2, numSeed);
    
    [pID,pN] = FDR(pValue,0.05)
    negLogPValue=zeros(numCommonState, numSeed);
    for i=1:numSeed
        for j=1:numCommonState
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
    end
end

% scatterplot all seeds/seedPairs and all states
if strcmp(measure,'TCPrctOverlapEachState')
    numSeed=6;
else
    numSeed=4;
end

numICC=length(find(ICCREML>-0.01))
indx=find(ICCREML>-0.01);

for j=1:numICC
    if j<=9
        figure(1)
        subplot(3,3,j)
        col=ceil(indx(j)/numCommonState)
        if rem(indx(j), numCommonState)~=0
            row=rem(indx(j), numCommonState)
        else
            row=numCommonState
        end
        scatter(data1(row,col,:), data2(row,col,:))
        
        if strcmp(measure,'stateFreq')
            xlim([-3 9]);
            ylim([-3 9]);
            set(gca,'xTick',-3:3:9, 'yTick', -3:3:9);
        elseif strcmp(measure,'TCPrctOverlapEachState')
            xlim([-0.5 1]);
            ylim([-0.5 1]);
            set(gca,'xTick',-0.5:0.5:1, 'yTick', -0.5:0.5:1);
       elseif strcmp(measure,'duration')
            xlim([-0.2 0.4]);
            ylim([-0.2 0.4]);
            set(gca,'xTick',-0.2:0.2:0.4, 'yTick', -0.2:0.2:0.4);
        else
            xlim([-0.5 1.5]);
            ylim([-0.5 1.5]);
            set(gca,'xTick',-0.5:0.5:1.5, 'yTick', -0.5:0.5:1.5);
        end
        
        if strcmp(measure,'TCPrctOverlapEachState')
            title(['seedPair', num2str(col), ' state', num2str(row)])
        else
            title(['seed', num2str(col), ' state', num2str(row)])
        end
        
        ylabel('Session 2')
        xlabel('Session 1')
        lsline
    elseif j<=18
        figure(2)
        subplot(3,3,j-9)
        col=ceil(indx(j)/numCommonState)
        if rem(indx(j), numCommonState)~=0
            row=rem(indx(j), numCommonState)
        else
            row=numCommonState
        end
        scatter(data1(row,col,:), data2(row,col,:))
        
        if strcmp(measure,'stateFreq')
            xlim([-3 9]);
            ylim([-3 9]);
            set(gca,'xTick',-3:3:9, 'yTick', -3:3:9);
        elseif strcmp(measure,'TCPrctOverlapEachState')
            xlim([-0.5 1]);
            ylim([-0.5 1]);
            set(gca,'xTick',-0.5:0.5:1, 'yTick', -0.5:0.5:1);
        elseif strcmp(measure,'duration')
            xlim([-0.2 0.4]);
            ylim([-0.2 0.4]);
            set(gca,'xTick',-0.2:0.2:0.4, 'yTick', -0.2:0.2:0.4);
        else
            xlim([-0.5 1.5]);
            ylim([-0.5 1.5]);
            set(gca,'xTick',-0.5:0.5:1.5, 'yTick', -0.5:0.5:1.5);
        end
        
        if strcmp(measure,'TCPrctOverlapEachState')
            title(['seedPair', num2str(col), ' state', num2str(row)])
        else
            title(['seed', num2str(col), ' state', num2str(row)])
        end
        
        ylabel('Session 2')
        xlabel('Session 1')
        lsline
    elseif j<=27
        figure(3)
        subplot(3,3,j-18)
        col=ceil(indx(j)/numCommonState)
        if rem(indx(j), numCommonState)~=0
            row=rem(indx(j), numCommonState)
        else
            row=numCommonState
        end
        scatter(data1(row,col,:), data2(row,col,:))
        
        if strcmp(measure,'stateFreq')
            xlim([-3 9]);
            ylim([-3 9]);
            set(gca,'xTick',-3:3:9, 'yTick', -3:3:9);
        elseif strcmp(measure,'TCPrctOverlapEachState')
            xlim([-0.5 1]);
            ylim([-0.5 1]);
            set(gca,'xTick',-0.5:0.5:1, 'yTick', -0.5:0.5:1);
        elseif strcmp(measure,'duration')
            xlim([-0.2 0.4]);
            ylim([-0.2 0.4]);
            set(gca,'xTick',-0.2:0.2:0.4, 'yTick', -0.2:0.2:0.4);
        else
            xlim([-0.5 1.5]);
            ylim([-0.5 1.5]);
            set(gca,'xTick',-0.5:0.5:1.5, 'yTick', -0.5:0.5:1.5);
        end
        
        if strcmp(measure,'TCPrctOverlapEachState')
            title(['seedPair', num2str(col), ' state', num2str(row)])
        else
            title(['seed', num2str(col), ' state', num2str(row)])
        end
        
        ylabel('Session 2')
        xlabel('Session 1')
        lsline
    else
        figure(4)
        subplot(3,3,j-27)
        col=ceil(indx(j)/numCommonState)
        if rem(indx(j), numCommonState)~=0
            row=rem(indx(j), numCommonState)
        else
            row=numCommonState
        end
        scatter(data1(row,col,:), data2(row,col,:))
        
        if strcmp(measure,'stateFreq')
            xlim([-3 9]);
            ylim([-3 9]);
            set(gca,'xTick',-3:3:9, 'yTick', -3:3:9);
        elseif strcmp(measure,'TCPrctOverlapEachState')
            xlim([-0.5 1]);
            ylim([-0.5 1]);
            set(gca,'xTick',-0.5:0.5:1, 'yTick', -0.5:0.5:1);
        elseif strcmp(measure,'duration')
            xlim([-0.2 0.4]);
            ylim([-0.2 0.4]);
            set(gca,'xTick',-0.2:0.2:0.4, 'yTick', -0.2:0.2:0.4);
        else
            xlim([-0.5 1.5]);
            ylim([-0.5 1.5]);
            set(gca,'xTick',-0.5:0.5:1.5, 'yTick', -0.5:0.5:1.5);
        end
        
        if strcmp(measure,'TCPrctOverlapEachState')
            title(['seedPair', num2str(col), ' state', num2str(row)])
        else
            title(['seed', num2str(col), ' state', num2str(row)])
        end
        
        ylabel('Session 2')
        xlabel('Session 1')
        lsline
    end
end
saveas(figure(1),[figDir,'scatter2SesAllSeedAllState_', measure, '(1).png'])
saveas(figure(2),[figDir,'scatter2SesAllSeedAllState_', measure, '(2).png'])
saveas(figure(3),[figDir,'scatter2SesAllSeedAllState_', measure, '(3).png'])
if strcmp(measure,'TCPrctOverlapEachState')
    saveas(figure(4),[figDir,'scatter2SesAllSeedAllState_', measure, '(4).png'])
end

% scatterplot ICC>0.4
if strcmp(measure,'TCPrctOverlapEachState')
    numSeed=6;
else
    numSeed=4;
end

numICC=length(find(ICCREML>0.4))
indx=find(ICCREML>0.4);
if numICC~=0
    figure(1)
    for j=1:numICC
        
        subplot(3,3,j)
        col=ceil(indx(j)/numCommonState)
        if rem(indx(j), numCommonState)~=0
            row=rem(indx(j), numCommonState)
        else
            row=numCommonState
        end
        scatter(data1(row,col,:), data2(row,col,:))
        
        if strcmp(measure,'stateFreq')
            xlim([-1 8]);
            ylim([-1 8]);
            set(gca,'xTick',-1:3:8, 'yTick', -1:3:8);
        else
            xlim([-0.5 1.5]);
            ylim([-0.5 1.5]);
            set(gca,'xTick',-0.5:0.5:1.5, 'yTick', -0.5:0.5:1.5);
        end
        
        if strcmp(measure,'TCPrctOverlapEachState')
            title(['seedPair', num2str(col), ' state', num2str(row)])
        else
            title(['seed', num2str(col), ' state', num2str(row)])
        end
        
        ylabel('Session 2')
        xlabel('Session 1')
        lsline
    end
    saveas(figure(1),[figDir,'scatter2Sessions_', measure, '.png'])
end


for k=1:numMeasure
    measure=char(measureList{k});
    figure(1)
    imagesc(ICCREML)
    colorbar
    title(['ICC',measure])
    ylabel('States')
    caxis([0 1])
    if strcmp(measure, 'TCPrctOverlapEachState')
        xlabel('Seed Pairs')
        set(gca,'xTick',1:6, 'yTick', 1:5);
        xlim([0.5 6.5]);
    else
        xlabel('Seeds')
        set(gca,'xTick',1:4, 'yTick', 1:5);
        xlim([0.5 4.5]);
    end
    ylim([0.5 5.5]);
    set(figure(1),'Colormap',mycmap)
    saveas(figure(1),[figDir,'ICC',measure,'.png'])
    
    figure(2)
    imagesc(mean1)
    colorbar
    title(['Session 1: ', measure])
    ylabel('States')
    if strcmp(measure, 'stateFreq')
        caxis([0 4])
    elseif strcmp(measure, 'subSeedPrct')
        caxis([0 1])
    elseif strcmp(measure, 'duration') || strcmp(measure, 'TCPrctOverlapEachState')
        caxis([0 0.2])
    end
    if strcmp(measure, 'TCPrctOverlapEachState')
        xlabel('Seed Pairs')
        set(gca,'xTick',1:6, 'yTick', 1:5);
        xlim([0.5 6.5]);
    else
        xlabel('Seeds')
        set(gca,'xTick',1:4, 'yTick', 1:5);
        xlim([0.5 4.5]);
    end
    ylim([0.5 5.5]);
    set(figure(2),'Colormap',mycmap)
    saveas(figure(2),[figDir,'Session1', measure,'.png'])
    
    figure(3)
    imagesc(mean2)
    colorbar
    title(['Session 2: ', measure])
    ylabel('States')
    if strcmp(measure, 'stateFreq')
        caxis([0 4])
    elseif strcmp(measure, 'subSeedPrct')
        caxis([0 1])
    elseif strcmp(measure, 'duration') || strcmp(measure, 'TCPrctOverlapEachState')
        caxis([0 0.2])
    end
    if strcmp(measure, 'TCPrctOverlapEachState')
        xlabel('Seed Pairs')
        set(gca,'xTick',1:6, 'yTick', 1:6);
        xlim([0.5 6.5]);
    else
        xlabel('Seeds')
        set(gca,'xTick',1:4, 'yTick', 1:6);
        xlim([0.5 4.5]);
    end
    ylim([0.5 6.5]);
    set(figure(3),'Colormap',mycmap)
    saveas(figure(3),[figDir,'Session2', measure,'.png'])
    
    
    figure(4)
    imagesc(negLogPValue)
    colorbar
    title(['Neglog10 of FDR corrected P Value: session dif on ', measure])
    
    ylabel('States')
    caxis([-5 5])
    if strcmp(measure, 'TCPrctOverlapEachState')
        xlabel('Seed Pairs')
        set(gca,'xTick',1:6, 'yTick', 1:5);
        xlim([0.5 6.5]);
    else
        xlabel('Seeds')
        set(gca,'xTick',1:4, 'yTick', 1:5);
        xlim([0.5 4.5]);
    end
    ylim([0.5 5.5]);
    saveas(figure(4),[figDir,'NegLogPValue', measure,'.png'])
    %close all
end




