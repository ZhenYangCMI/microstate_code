


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
figDir='/Users/zhenyang/Desktop/Zhen/figs/5_3_13/';

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

% scatterplot the two transitions with ICC>0.4

close all

for i=1:numSeed
numICC=length(find(ICCREML(:,:,i)>0.4))
    indx=find(ICCREML(:,:,i)>0.4);
    figure(i)
    for j=1:numICC
        subplot(2,3,j)
        col=ceil(indx(j)/numState)
        if rem(indx(j), numState)~=0
        row=rem(indx(j), numState)
        else
            row=numState
        end
        scatter(data1(row,col,:,i), data2(row,col,:,i))
        xlim([-0.1 0.3]);
        ylim([-0.1 0.3]);
        set(gca,'xTick',-0.1:0.1:0.3, 'yTick', -0.1:0.1:0.3);
        ylabel('Session 2')
        xlabel('Session 1')
        title(['Seed', num2str(i), ': state', num2str(row),' (t) to state', num2str(col),' (t+1)'])
        lsline
    end
    saveas(figure(i),[figDir,'scatter2Sessions_HighICC_TM_seed_replot',num2str(i),'.png'])
end



% assign 1 to the diagonal for plotting purpose
for i=1:numSeed
    dataSes1=squeeze(data5(:,:,i));
    dataSes2=squeeze(data6(:,:,i));
    % remove the probability of entering or exiting its own state
    for m=1:numState
        dataSes1(m,m)=1;
    end
    for n=1:numState
        dataSes2(n,n)=1;
    end
    data5(:,:,i)=dataSes1;
    data6(:,:,i)=dataSes2;
end



[pID,pN] = FDR(pValue,0.05)
negLogPValue=zeros(numState, numSeed);
for i=1:numSeed
    for j=1:numState
        for k=1:numState
            if pValue(j,k,i) <= pID
                if mean1(j,k,i)-mean2(j,k,i)>0;
                    negLogPValue(j,k,i)=(-1)*log10(pValue(j,k,i));
                else
                    negLogPValue(j,k,i)=log10(pValue(j,k,i));
                end
            else
                negLogPValue(j,k,i)=0;
            end
        end
    end
end
negLogPValue

for i=1:numSeed
    dataPlot3=data3(:,:,i);
    dataPlot4=data4(:,:,i);
    
    dataPlot5=data5(:,:,i);
    dataPlot6=data6(:,:,i);
    
    figure(1)
    subplot(2,2,i)
    imagesc(dataPlot3)
    h=colorbar;
    set(h, 'YScale', 'log', 'ytick',[1e-3,1e-2,1e-1,1e0],  'YMinorTick', 'off')
    my_scale = [1e-3 1e0];
    caxis(my_scale)
    set(gca,'xTick',1:5, 'yTick', 1:5);
    xlim([0.5 5.5]);
    ylim([0.5 5.5]);
    set(figure(1),'Colormap',mycmap3)
    ylabel('States at time t')
    xlabel('States at time t+1')
    title(['Transition Matrix of Markov Chain Session1 seed', num2str(i)])
    saveas(figure(1),[figDir,'Session1', measure1, '.png'])
    
    figure(2)
    subplot(2,2,i)
    imagesc(dataPlot4)
    h=colorbar;
    set(h, 'YScale', 'log', 'ytick',[1e-3,1e-2,1e-1,1e0],  'YMinorTick', 'off')
    my_scale = [1e-3 1e0];
    caxis(my_scale)
    set(gca,'xTick',1:5, 'yTick', 1:5);
    xlim([0.5 5.5]);
    ylim([0.5 5.5]);
    set(figure(2),'Colormap',mycmap3)
    ylabel('States at time t')
    xlabel('States at time t+1')
    title(['Transition Matrix of Markov Chain Session2 seed', num2str(i)])
    saveas(figure(2),[figDir,'Session2', measure1,'.png'])
    
    figure(3)
    subplot(2,2,i)
    imagesc(ICCREML(:,:,i))
    colorbar
    title(['ICC',measure, ' seed', num2str(i)])
    ylabel('States at time t')
    xlabel('States at time t+1')
    set(gca,'xTick',1:5, 'yTick', 1:5);
    caxis([0 1])
    xlim([0.5 5.5]);
    ylim([0.5 5.5]);
    set(figure(3),'Colormap',mycmap)
    saveas(figure(3),[figDir,'ICC',measure,'.png'])
    
    figure(4)
    subplot(2,2,i)
    imagesc(negLogPValue(:,:,i))
    colorbar
    title(['Neglog10: session dif on ', measure, ' seed', num2str(i), '(FDR corrected)'])
    ylabel('States at time t')
    xlabel('States at time t+1')
    set(gca,'xTick',1:5, 'yTick', 1:5);
    xlim([0.5 5.5]);
    ylim([0.5 5.5]);
    caxis([-5 5])
    saveas(figure(4),[figDir,'NegLogPValue', measure,'.png'])
    
    figure(5)
    subplot(2,2,i)
    imagesc(dataPlot5)
    h=colorbar;
    set(h, 'YScale', 'log', 'ytick',[1e-3,1e-2,1e-1,1e0],  'YMinorTick', 'off')
    my_scale = [1e-3 1e0];
    caxis(my_scale)
    set(gca,'xTick',1:5, 'yTick', 1:5);
    xlim([0.5 5.5]);
    ylim([0.5 5.5]);
    set(figure(5),'Colormap',mycmap3)
    ylabel('States at time t')
    xlabel('States at time t+1')
    title(['TM to dif state Session1 seed', num2str(i)])
    saveas(figure(5),[figDir,'Session1', measure2, '.png'])
    
    figure(6)
    subplot(2,2,i)
    imagesc(dataPlot6)
    h=colorbar;
    set(h, 'YScale', 'log', 'ytick',[1e-3,1e-2,1e-1,1e0],  'YMinorTick', 'off')
    my_scale = [1e-3 1e0];
    caxis(my_scale)
    set(gca,'xTick',1:5, 'yTick', 1:5);
    xlim([0.5 5.5]);
    ylim([0.5 5.5]);
    set(figure(6),'Colormap',mycmap3)
    ylabel('States at time t')
    xlabel('States at time t+1')
    title(['TM to dif state Session2 seed', num2str(i)])
    saveas(figure(6),[figDir,'Session2', measure2, '.png'])
end



