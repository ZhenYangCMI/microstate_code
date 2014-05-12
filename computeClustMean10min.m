
clear
clc
close all

sessionList={'session1'};
numSession=length(sessionList)
dataLength='all_10min';
winSizeList=[136];

for k=1:length(winSizeList)
    winSize=winSizeList(k);
    
    for i=1:numSession
        session=char(sessionList{i});
        resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength, filesep, session,filesep];
                
        tmp=load([resultDir,'/zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'win',num2str(winSize),'.mat']);
        winAllSubAllSeedFullCor=tmp.zWinFullCorLasso;
        numROI=156;
        
        indx=load([resultDir, 'clustIndxWinAllSeeds_FullCorLasso_',session,'_10min_win', num2str(winSize),'.txt']);
        numClust=length(unique(indx));
        text=sprintf('Working on surrogate%d %s_%dclusters_win',k,session, numClust);
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
            save([resultDir, sprintf('clustMean/clusterMean_%dclusters_%s_normWin_win%d.mat',numClust,session,winSize)],'finalMeanWinOfClust')
            
            
        end
    end
end
