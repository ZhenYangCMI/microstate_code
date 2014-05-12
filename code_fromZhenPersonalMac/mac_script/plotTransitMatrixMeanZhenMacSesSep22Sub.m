


clear
clc
close all


numSeed=4;
numSub=22;
numState=5;
session='session1'

dataLength='all_10min';

colors=[0,100/255,0; 34/255,139/255,34/255;50/255,205/255,50/255;152/255,251/255,152/255;240/255,255/255,240/255;1,1,1];
colorMap = colorRamp(colors, 32)

resultDir=['/Users/zhenyang/Desktop/Zhen/results/all_10min/',session,filesep];
figDir='/Users/zhenyang/Desktop/Zhen/figs/5_24_13/';

measure='probTM'
tmp1=load([resultDir, measure,'_',session,'_recode.mat']);
data1=tmp1.(measure);

measure1='meanProbTM';
tmp2=load([resultDir, measure1,'_',session,'_recode.mat']);
data2=tmp2.(measure1);

measure2='meanProbTMDifState';
tmp3=load([resultDir, measure2,'_',session,'_recode.mat']);
data3=tmp3.(measure2);




% assign 1 to the diagonal for plotting purpose
for i=1:numSeed
    tmp=squeeze(data3(:,:,i));
    
    % remove the probability of entering or exiting its own state
    for m=1:numState
        tmp(m,m)=1;
    end
       data3(:,:,i)=tmp;
    end


for i=1:numSeed
        figure(i)
    imagesc(data3(:,:,i))
   h=colorbar;
    set(h, 'YScale', 'log', 'ytick',[1e-3,1e-2,1e-1,1e0],  'YMinorTick', 'off')
       my_scale = [1e-3 1e0];
    caxis(my_scale)
        colorbar('off')
    axis square
    set(gca,'xTick',0.5:5.5, 'yTick', 0.5:5.5,'xTickLabel',[],'yTickLabel',[]);
   
        grid on
    set(gca,'GridLineStyle','-')
    set(gca,'LineWidth',2)
    xlim([0.5 5.5]);
    ylim([0.5 5.5]);
    set(figure(i),'Colormap',colorMap)
    %saveas(figure(i),[figDir,measure,'_seed',num2str(i),'_session1.png'])
end



