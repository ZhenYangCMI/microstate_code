% This script compute the mean FD for eash state each seed
clear
clc

session='session1';
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;

numSeed=4;
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
%finalmeanWinFDPower=repmat(meanWinFDPower, 4, 1);

resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/', session,'/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/', session,'/motion/'];

indx=load([resultDir, 'clustIndxNormWinAllSeeds_FullCorLasso_',session,'_10min.txt']);

numWinPerSeed=5984;

for k=1:numSeed
    indxSeed=indx((1+numWinPerSeed*(k-1)):numWinPerSeed*k);
    numClust=length(unique(indxSeed));
    clustMeanFD=zeros(numWinPerSeed, numClust);
    numWinInClust=zeros(numClust);
    group=zeros(numClust);
    for i=1:numClust
        indxList=find(indxSeed==i);
        numWinInClust(i)=length(indxList);
        
        tmp=meanWinFDPower(indxList);
        for j=1:numWinInClust(i)
            clustMeanFD(j, i)=tmp(j,1);
        end
    end
    clustMeanFD(find(clustMeanFD==0))=[];
    clustMeanFD=clustMeanFD';
    numWinInClust1=numWinInClust(1);
    numWinInClust2=numWinInClust(2);
    numWinInClust3=numWinInClust(3);
    numWinInClust4=numWinInClust(4);
    numWinInClust5=numWinInClust(5);
    
    group1=repmat(1, numWinInClust1,1);
    group2=repmat(2, numWinInClust2,1);
    group3=repmat(3, numWinInClust3,1);
    group4=repmat(4, numWinInClust4,1);
    group5=repmat(5, numWinInClust5,1);
        group=[group1; group2; group3; group4; group5];
    
    %p = vartestn(clustMeanFD,group,'TestType','LeveneAbsolute')
    
    close all
    anova1(clustMeanFD, group)
    kruskalwallis(clustMeanFD, group)
    saveas(figure(1), [figDir, 'seed', num2str(k), '_anovaTable.png'])
    saveas(figure(2), [figDir, 'seed', num2str(k), '_boxplot.png'])
    
end

% clust1MeanFD=clustMeanFD(1:5620, 1);
% clust2MeanFD=clustMeanFD(5621:5621+5602-1,1);
% clust3MeanFD=clustMeanFD(5621+5602:5621+5602+4822-1, 1);
% clust4MeanFD=clustMeanFD(5621+5602+4822:5621+5602+4822+4632-1, 1);
% clust5MeanFD=clustMeanFD(5621+5602+4822+4632:5621+5602+4822+4632+3260-1,1);
%
% [h, p, ci, stats]=ttest2(clust3MeanFD, clust5MeanFD, 0.05,'both', 'unequal')
% [h, p, ci, stats]=ttest2(clust1MeanFD, clust2MeanFD, 0.05,'both', 'unequal')
%
%
% save([resultDir, 'FDPower/clustMeanFD.txt'], 'clustMeanFD','-ascii', '-double', '-tabs')
