
clear
clc

resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/session1/ses1_ward_euclidean/'];

for i=2:2:16
    numClust=i;
    tmp=load([resultDir,'indx_',num2str(numClust),'clusters.mat']);
    indx=tmp.indx;
    numWinInClust=countWinInClust(indx)
    save([resultDir,'numWinInClust_',num2str(numClust),'clusters.mat'], 'numWinInClust')
end 



