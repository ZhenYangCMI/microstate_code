clear
clc
close all

resultDir='/Users/zhenyang/Desktop/Zhen/results/all_10min/2sesConcate/';
figDir='/Users/zhenyang/Desktop/Zhen/figs/5_22_13/';
numSub=21;
measure='duration';
seed=4
state=4

% scatterplot the duration_seed4_state4 and DKEFS measure
numSeed=4;
numSubExclude=5;
numState=5;
% sub 1-4 and 11 were removed, after removing each one the order will be different, so the Indx is different from the orignal order
measure='duration';
subIDList=[1,1,1,1,7];
NPList={'DKEFS_DESIGN_SWITCHING_SCALED', 'DKEFS_VERB_CATEGORYSWITCH_SCALED', 'DKEFS_TOTALWEIGHTACHIEVESCORE_SCALED',...
    'DKEFS_COLORWORD_INHIBITION_SCALED','WASI_PERF','WASI_VERB','WASI_FULL_4'};
NPdata1=load([resultDir,'DKEFS_WASI_16sub_2sesCmb_session1.txt']);
NPdata2=load([resultDir,'DKEFS_WASI_16sub_2sesCmb_session2.txt']);

tmp1=load([resultDir, measure,'_session1.mat']);
data1=tmp1.(measure);
tmp2=load([resultDir, measure,'_session2.mat']);
data2=tmp2.(measure);
data16sub1=data1;
data16sub2=data2;
for p=1:numSubExclude
        subID=subIDList(p);
    data16sub1(:,:,subID)=[];
    data16sub2(:,:,subID,:)=[];
end

lineWidth=2.5;
duration1=squeeze(data16sub1(state,seed,:));
duration2=squeeze(data16sub2(state,seed,:));
NP1=squeeze(NPdata1(:,3));
NP2=squeeze(NPdata2(:,3));
fcolor=[0,0,255];
fcolor=fcolor/255
ecolor=[0,0,1];
s=300
figure(2)
scatter(duration1, NP1, s,'o','MarkerEdgeColor',ecolor,'MarkerFaceColor',fcolor,'LineWidth',1)
%set(gca,'xTick',-0.1:0.1:0.3, 'yTick', 2:4:18, 'XTickLabel',[],'YTickLabel',[]);
% ylim dosen't work, always set back to ~-5 to 21. mannually edited it to 0
% to 20, mannually edited the linewidth of axis to 1.2
set(gca,'xlim',[-0.1 0.3],'ylim',[0 20],'xTick',-0.1:0.1:0.3, 'yTick', 0:5:20, 'XTickLabel',[],'YTickLabel',[]);
h=lsline;
set(h,'LineWidth',lineWidth)
set(gca,'LineWidth',lineWidth)
axis square
saveas(figure(2),[figDir, 'scatterPhenotypeDyn_session1_1.png'])

figure(3)

scatter(duration2, NP2, s,'o','MarkerEdgeColor',ecolor,'MarkerFaceColor',fcolor,'LineWidth',1)
%set(gca,'xTick',-0.1:0.1:0.3, 'yTick', 2:4:18, 'XTickLabel',[],'YTickLabel',[]);
% mannually edited the linewidth of axis to 1.5
set(gca,'xlim',[-0.1 0.3],'ylim',[0 20],'xTick',-0.1:0.1:0.3, 'yTick', 0:5:20, 'XTickLabel',[],'YTickLabel',[]);
h=lsline;
set(h,'LineWidth',lineWidth)
set(gca,'LineWidth',lineWidth)
axis square
saveas(figure(3),[figDir, 'scatterPhenotypeDyn_session2_1.png'])



% Plot scatter plot of duration_seed4_state4 between 2 sessions
tmp1=load([resultDir, measure,'_session1.mat']);
data1=tmp1.(measure);
tmp2=load([resultDir, measure,'_session2.mat']);
data2=tmp2.(measure);

%fcolor=[238,130,238];
fcolor=[255,0,0];
fcolor=fcolor/255
ecolor=[1,0,0];
s=300
lineWidth=2.5;
d1=squeeze(data1(state,seed,:));
d2=squeeze(data2(state,seed,:));

figure(1)
scatter(d1, d2, s, 'o','MarkerEdgeColor',ecolor,'MarkerFaceColor',fcolor,'LineWidth',1)
xlim([-0.1 0.3]);
ylim([-0.1 0.3]);
set(gca,'xTick',-0.1:0.1:0.3, 'yTick', -0.1:0.1:0.3, 'XTickLabel',[],'YTickLabel',[]);
h=lsline;
set(h,'LineWidth',lineWidth)
set(gca,'LineWidth',lineWidth)
axis square
saveas(figure(1),[figDir, 'scatterSes1Ses2_duration_2sesCmb21Sub1.png'])