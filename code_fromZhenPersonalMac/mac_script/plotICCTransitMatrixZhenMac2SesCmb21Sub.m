


clear
clc
close all


numSeed=4;
numSub=21;
numState=5;


dataLength='all_10min';
load('MyColormaps','mycmap') % [0 1]
load('MyColormaps2','mycmap2') % > 1
load('MyColormaps3','mycmap3') % log scale, white largest value

resultDir='/Users/zhenyang/Desktop/Zhen/results/all_10min/2sesConcate/';
figDir='/Users/zhenyang/Desktop/Zhen/figs/5_24_13/';

measure='probTM'
tmp1=load([resultDir, measure,'_session1.mat']);
data1=tmp1.(measure);
tmp2=load([resultDir, measure,'_session2.mat']);
data2=tmp2.(measure);

measure1='meanProbTM';
tmp3=load([resultDir, measure1,'_session1.mat']);
data3=tmp3.(measure1);
tmp4=load([resultDir, measure1,'_session2.mat']);
data4=tmp4.(measure1);

measure2='meanProbTMDifState';
tmp5=load([resultDir, measure2,'_session1.mat']);
data5=tmp5.(measure2);
tmp6=load([resultDir, measure2,'_session2.mat']);
data6=tmp6.(measure2);


ICCREML=zeros(numState,numState, numSeed);
ICCIPN=zeros(numState,numState,numSeed);
for i=1:numSeed
    disp(['Work on seed ', num2str(i)])
    for j=1:numState
        for k=1:numState
            session1Data=squeeze(data1(j,k,:,i));
            session2Data=squeeze(data2(j,k,:,i));
            data=[session1Data,session2Data];
            Y=[session1Data;session2Data];
            time = [ones(numSub,1);2*ones(numSub,1)];
            sID=[[1:numSub]';[1:numSub]'];
            [ ICC1, idx_fail] = do_ICC(Y, time, [], [], sID);
            ICCREML(j,k,i)=ICC1;
            ICC2 = IPN_icc(data,1,'single');
            ICCIPN(j,k,i)=ICC2;
            [h,p]=ttest(session1Data, session2Data);
            pValue(j,k,i)=p;
        end
    end
end
ICCREML
ICCIPN
pValue


% assign 1 to the diagonal for plotting purpose
for i=1:numSeed
    tmp=squeeze(ICCREML(:,:,i));
    
    % remove the probability of entering or exiting its own state
    for m=1:numState
        tmp(m,m)=1;
    end
       ICCREML(:,:,i)=tmp;
    end

colors=[1,0.5,0;1,0.7,0;1,1,0; 255/255,250/255,205/255;255/255,255/255,224/255;1,1,1];
colorMap = colorRamp(colors, 32)

for i=1:numSeed
        figure(i)
    imagesc(ICCREML(:,:,i))
        caxis([0 1])
        %hcb=colorbar;
%     set(hcb,'YTick',[0,0.5,1])
    axis square
    set(gca,'xTick',0.5:5.5, 'yTick', 0.5:5.5,'xTickLabel',[],'yTickLabel',[]);
        grid on
    set(gca,'GridLineStyle','-')
    set(gca,'LineWidth',2)
    xlim([0.5 5.5]);
    ylim([0.5 5.5]);
    set(figure(i),'Colormap',colorMap)
    %saveas(figure(i),[figDir,'ICC',measure,'_seed',num2str(i),'.png'])
end



