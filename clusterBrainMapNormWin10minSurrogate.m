
clear
clc
close all

sessionList={'session1'};
numSession=length(sessionList)
dataLength='all_10min';
numSurrogate=100;

for k=1:numSurrogate
    
    for i=1:numSession
        session=char(sessionList{i});
        resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength, filesep, session,filesep];
        figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep, session,'/surrogateData/'];
        
        tmp=load([resultDir,'surrogate/zNormWindows/zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'_surrogate_',num2str(k),'.mat']);
        winAllSubAllSeedFullCor=tmp.zWinFullCorLasso;
        numROI=156;
        
        indx=load([resultDir, 'clustIndxNormWinAllSeeds_FullCorLasso_',session,'_10min.txt']);
        numClust=length(unique(indx));
        text=sprintf('Working on surrogate%d %s_%dclusters',k,session, numClust);
        disp(text)
        
        for i=1:numClust
            indxList=find(indx==i);
            numWinInClust=length(indxList);
            allWinInClust=zeros(numWinInClust,numROI);
            for j=1:numWinInClust
                for m=1:numROI
                    allWinInClust(j,m)=winAllSubAllSeedFullCor(indxList(j), m);
                end
            end
            meanWinOfClust=mean(allWinInClust);
            for m=1:numROI
                finalMeanWinOfClust(i,m)=meanWinOfClust(1,m);
            end
            disp('Cluster centroids computed. Start saving file.')
            save([resultDir, sprintf('surrogate/clustMean/clusterMean_%dclusters_%s_normWin_surrogate_%d.mat',numClust,session,k)],'finalMeanWinOfClust')
        end
    end
end
