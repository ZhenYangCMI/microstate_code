clear
clc

session='session1';
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;
covType='GSR'

winWidth=69; % in TR
step=3; % in TR
numVol=884;
numWin=floor((numVol-winWidth)/step)+1;
numSub=length(subList);
numWinAllSub=numWin*numSub;
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];

winFDPower=zeros(winWidth,numWinAllSub);
t=0
for k=1:numSub
    sub=subList{k};
    disp (['Working on sub ', char(sub),' ......'])
    subDir=[analyDir,'data/645/all_10min/preprocessed/RealignParameter/',char(sub)];
    FDPower=load([subDir, '/FD_Power_', char(sub), '.txt']);
    
    for i=1:step:(numVol-winWidth+1)
        t=t+1
        winFDPower(:, t)=FDPower(i:i+winWidth-1, 1);
    end
end

meanWinFDPower=mean(winFDPower)';
finalmeanWinFDPower=repmat(meanWinFDPower, 4, 1);

resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/', covType, '/', session,'/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/', covType, '/', session,'/motion/'];

indx=load([resultDir, 'clustIndxNormWinAllSeeds_FullCorLasso_',session,'_10min.txt']);
numClust=length(unique(indx));
clustMeanFD=zeros(5620, numClust);
for i=1:numClust
    indxList=find(indx==i);
    numWinInClust=length(indxList)
    tmp=finalmeanWinFDPower(indxList);
    for j=1:numWinInClust
        clustMeanFD(j, i)=tmp(j,1);
    end
end

for k=1:numClust
    a=clustMeanFD(:, k);
    a(find(a==0))=[];
    c=length(find(a>0.2))/length(a)*100;
    b(k, 1)=median(a)
    percent(k, 1)=c;
end

lineWidth=2
figure(1)
h=bar(percent,0.4) % TCPrctOverlapStateSum used 0.6, other 2 measures used 0.4 to make them the same width
set(gca,'LineWidth',lineWidth)
set(h,'LineWidth',lineWidth,'Facecolor',[0,0,1],'EdgeColor',[0 0 1])
set(gca,'xTick',1:5, 'yTick', 0:10:40, 'box','off', 'xTickLabel',[],'yTickLabel',[]);
ylim([0 40])
saveas(figure(1), sprintf('/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/GSR/%s/percentWinScrub_%s.png', session, session))

% The following is for anova analysis
%clustMeanFD(find(clustMeanFD==0))=[];
%clustMeanFD=clustMeanFD';

clustMeanFD(find(clustMeanFD==0))=NaN;

close all
lineWidth=2.5
figure(2)
h=boxplot(clustMeanFD, 'whisker', 3, 'outliersize', 5)
    %http://alex.bikfalvi.com/research/advanced_matlab_boxplot/
        set(gca,'LineWidth',2)
        set(h, 'LineWidth', lineWidth)
     ylim([-0.2 0.8]);
        set(gca, 'yTick', -0.2:0.2:0.8, 'xTick', 1:5, 'xTickLabel',[], 'yTickLabel',[], 'box','off');
        set(gca, 'xTickLabel',[], 'yTickLabel',[], 'box','off')
        saveas(figure(2), sprintf('/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/GSR/%s/boxplotMotion_%s.png', session, session))




group1=repmat(1, 5620, 1);
group2=repmat(2, 5602, 1);
group3=repmat(3, 4822, 1);
group4=repmat(4, 4632, 1);
group5=repmat(5, 3260, 1);
group=[group1; group2; group3; group4; group5];

%p = vartestn(clustMeanFD,group,'TestType','LeveneAbsolute')


close all
anova1(clustMeanFD, group)
kruskalwallis(clustMeanFD, group)
% saveas(figure(1), [figDir, 'anovaTable.png'])
% saveas(figure(2), [figDir, 'boxplot.png'])


clust1MeanFD=clustMeanFD(1:5620, 1);
clust2MeanFD=clustMeanFD(5621:5621+5602-1,1);
clust3MeanFD=clustMeanFD(5621+5602:5621+5602+4822-1, 1);
clust4MeanFD=clustMeanFD(5621+5602+4822:5621+5602+4822+4632-1, 1);
clust5MeanFD=clustMeanFD(5621+5602+4822+4632:5621+5602+4822+4632+3260-1,1);

[h, p, ci, stats]=ttest2(clust3MeanFD, clust5MeanFD, 0.05,'both', 'unequal')
[h, p, ci, stats]=ttest2(clust1MeanFD, clust2MeanFD, 0.05,'both', 'unequal')


save([resultDir, 'FDPower/clustMeanFD.txt'], 'clustMeanFD','-ascii', '-double', '-tabs')
