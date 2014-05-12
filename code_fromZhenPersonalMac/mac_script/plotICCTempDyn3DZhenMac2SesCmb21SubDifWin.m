clear
clc
close all


numSub=21;
numState=4;
winSize=136;
numSub=21;

dataLength='all_10min';

resultDir='/Users/zhenyang/Desktop/Zhen/results/all_10min/2sesConcate/';
figDir='/Users/zhenyang/Desktop/Zhen/figs/7_1_13/';

%measureList={'subSeedPrct', 'stateFreq', 'duration', 'TCPrctOverlapEachState'};
measureList={'duration'}
numMeasure=length(measureList);

colors=[1,0.5,0;1,0.7,0;1,1,0; 255/255,250/255,205/255;255/255,255/255,224/255;1,1,1];
colorMap = colorRamp(colors, 32)

for k=1:numMeasure
    measure=char(measureList{k});
    if strcmp(measure,'TCPrctOverlapEachState')
        numSeed=6;
    else
        numSeed=4;
    end
    % compute ICC and perform ttest on % time spent in a state
    tmp1=load([resultDir, measure,'_session1_win', num2str(winSize), '_', num2str(numSub), 'sub.mat']);
    data1=tmp1.(measure);
    tmp2=load([resultDir, measure,'_session2_win', num2str(winSize), '_', num2str(numSub), 'sub.mat']);;
    data2=tmp2.(measure);
    
    ICCREML=zeros(numState,numSeed);
    ICCIPN=zeros(numState,numSeed);
    for i=1:numSeed
        disp(['Work on seed ', num2str(i)])
        for j=1:numState
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
    mean1=reshape(mean1Reshap, numState, numSeed);
    data2Reshape=reshape(data2, [], numSub)';
    mean2Reshap=mean(data2Reshape);
    mean2=reshape(mean2Reshap, numState, numSeed);
    
    [pID,pN] = FDR(pValue,0.05)
    negLogPValue=zeros(numState, numSeed);
    for i=1:numSeed
        for j=1:numState
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
    
    figure(1)
    imagesc(ICCREML)
    caxis([0 1])
    %     hcb=colorbar;
    %     set(hcb,'YTick',[0,0.5,1])
    %     title(['ICC',measure])
    %     ylabel('States')
    
    grid on
    set(gca,'GridLineStyle','-')
    set(gca,'LineWidth',2)
    if strcmp(measure, 'TCPrctOverlapEachState')
        set(gca,'xTick',0.5:6.5, 'yTick', 0.5:5.5,'xTickLabel',[],'yTickLabel',[]);
        xlim([0.5 6.5]);
    else
        set(gca,'xTick',0.5:4.5, 'yTick', 0.5:4.5, 'xTickLabel',[],'yTickLabel',[]);
        xlim([0.5 4.5]);
    end
    ylim([0.5 4.5]);
    set(figure(1),'Colormap',colorMap)
    saveas(figure(1),[figDir,sprintf('ICC_%s_win%d_%dsub.png',measure,winSize,numSub)])
    
    figure(2)
    scatter(data1(4,4, :), data2(4,4,:))
    lsline
    saveas(figure(2),[figDir,sprintf('scatter2Ses_%s_win%d_%dsub.png',measure,winSize,numSub)])
end
