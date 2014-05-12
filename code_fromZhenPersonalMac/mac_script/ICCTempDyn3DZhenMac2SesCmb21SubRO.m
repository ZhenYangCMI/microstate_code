clear
clc
close all

numState=5;

dataLength='all_10min';
load('MyColormaps2','mycmap2')
load('MyColormaps','mycmap')


resultDir='/Users/zhenyang/Desktop/Zhen/results/all_10min/2sesConcate/';
figDir='/Users/zhenyang/Desktop/Zhen/figs/5_3_13/';


measure='duration';
numSeed=4;
% compute ICC and perform ttest on % time spent in a state
tmp1=load([resultDir, measure,'_session1.mat']);
data1=tmp1.(measure);
tmp2=load([resultDir, measure,'_session2.mat']);
data2=tmp2.(measure);

% seed4State3, remove the 2 outliers, recompute the ICC and
% rescatterplot
ICCREML=zeros(numState,numSeed);
numSub=19;
for i=4 %seed num
    disp(['Work on seed ', num2str(i)])
    for j=3 % state num
        session1Data=squeeze(data1(j,i,:));
        session2Data=squeeze(data2(j,i,:));
        session1Data(12)=[];
        session1Data(14)=[];
        session2Data(12)=[];
        session2Data(14)=[];
        Y=[session1Data;session2Data];
        time = [ones(numSub,1);2*ones(numSub,1)];
        sID=[[1:numSub]';[1:numSub]'];
        [ ICC1, idx_fail] = do_ICC(Y, time, [], [], sID);
        ICCREML(j,i)=ICC1;
        
        ICCREML
        
        figure(1)
        scatter(session1Data, session2Data, 100)
        xlim([-0.2 0.6]);
        ylim([-0.2 0.6]);
        set(gca,'xTick',-0.2:0.2:0.6, 'yTick', -0.2:0.2:0.6);
        h=lsline
        set(h,'LineWidth', 1.5)
        ylabel('Session 2')
        xlabel('Session 1')
        saveas(figure(1),[figDir,sprintf('scatter2seed%dstate%d_wo.png', i, j)])
    end
end

% seed2State2, remove the four outliers: sub 16 **, recompute the ICC and
% rescatterplot
ICCREML=zeros(numState,numSeed);
numSub=17;
for i=2 %seed num
    disp(['Work on seed ', num2str(i)])
    for j=2 % state num
        session1Data=squeeze(data1(j,i,:));
        session2Data=squeeze(data2(j,i,:));
        session1Data(find(session2Data==0|session2Data>0.2))=[];
        session2Data(find(session2Data==0|session2Data>0.2))=[];
        Y=[session1Data;session2Data];
        time = [ones(numSub,1);2*ones(numSub,1)];
        sID=[[1:numSub]';[1:numSub]'];
        [ ICC1, idx_fail] = do_ICC(Y, time, [], [], sID);
        ICCREML(j,i)=ICC1;
        
        ICCREML
        
        figure(1)
        scatter(session1Data, session2Data, 100)
        xlim([-0.2 0.4]);
        ylim([-0.2 0.4]);
        set(gca,'xTick',-0.2:0.2:0.4, 'yTick', -0.2:0.2:0.4);
        h=lsline
        set(h,'LineWidth', 1.5)
        ylabel('Session 2')
        xlabel('Session 1')
        saveas(figure(1),[figDir,sprintf('scatter2seed%dstate%d_wo.png', i, j)])
    end
end

% seed4State4, remove the outlier sub 16, recompute the ICC and
% rescatterplot
ICCREML=zeros(numState,numSeed);
numSub=20;
for i=4 %seed num
    disp(['Work on seed ', num2str(i)])
    for j=4 % state num
        session1Data=squeeze(data1(j,i,:));
        session2Data=squeeze(data2(j,i,:));
        session1Data(16)=[];
        session2Data(16)=[];
        Y=[session1Data;session2Data];
        time = [ones(numSub,1);2*ones(numSub,1)];
        sID=[[1:numSub]';[1:numSub]'];
        [ ICC1, idx_fail] = do_ICC(Y, time, [], [], sID);
        ICCREML(j,i)=ICC1;
        
        ICCREML
        
        figure(1)
        scatter(session1Data, session2Data, 100)
        xlim([-0.2 0.4]);
        ylim([-0.2 0.4]);
        set(gca,'xTick',-0.2:0.2:0.4, 'yTick', -0.2:0.2:0.4);
        h=lsline
        set(h,'LineWidth', 1.5)
        ylabel('Session 2')
        xlabel('Session 1')
        saveas(figure(1),[figDir,sprintf('scatter2seed%dstate%d_wo.png', i, j)])
    end
end