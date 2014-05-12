function [ R2 PseudoF ] = clustNumEstimate( numClust session )
%UNTITLED Summary of this function goes here
%   input: 1. session: string, e.g. session='session1'

resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',session,'/'];

tmp1=load([resultDir,'winAllSubAllSeed_partialCor_645_',session,'.mat']);
winAllSubAllSeedPartialCor=tmp1.partWinAllSubAllSeed;

tmp2=load([resultDir,session,'_ward_euclidean/computeClustIndx/indxList_2_to_100_clusters.mat']);
indxList=tmp2.indxList;
disp ('The window and index files are loaded successfully.')

numWin=size(winAllSubAllSeedPartialCor, 1);
numROI=size(winAllSubAllSeedPartialCor, 2);

% compute center of all windows
meanWinAllSubAllSeedPartialCor=mean(winAllSubAllSeedPartialCor);

% compute the sum distance of each window to the ground mean window
T=0;
for i=1:numWin
    for j=1:numROI
        T=T+(winAllSubAllSeedPartialCor(i,j)-meanWinAllSubAllSeedPartialCor(j))^2;
    end
end

% compute the PG: sum squares within each cluster

for m=2:100
    indx=indxList(:,m-1);
    numClust=m;
    for n=1:numClust
        clustIndx=find(indx==n);
        numWinInClust=length(clustIndx);
        allWinInClust=zeros(numWinInClust,numROI);
        for p=1:numWinInClust
            for q=1:numROI
                allWinInClust(p,q)=winAllSubAllSeedPartialCor(clustIndx(p), q);
            end
        end
        meanWinOfClust=mean(allWinInClust);
        k=0;PG=zeros(numWin*numROI);
        for p=1:numWin
            for q=1:numROI
                k=k+1;
                PG(k)=(allWinInClust(p,q)-meanWinOfClust(p,q))^2;
            end
        end
        sumPG=sum(PG);
        PGAllClust(numClust)=
        for q=1:numROI
            finalMeanWinOfClust(n,q)=meanWinOfClust(1,q);
        end
        
        %R^2 method
        R2(m-1)=1-PG(m-1)/T;
        
        %pseudo-F statistic method
        pseudoF(m-1)=((T-PG(m-1))/(m-1))/(PG(m-1)/(numWin-m));
        
        clear
    end
    
