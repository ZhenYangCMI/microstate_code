
clear
clc
close all

sessionList={'session1','session2'};
numSession=length(sessionList)
dataLength='all_10min';
for i=1:numSession
    session=char(sessionList{i});
    resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength, filesep, session,'/'];
    maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/'];
    figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep, session,'/corBrainMapAvg/'];
    
    tmp=load([resultDir,'zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'.mat']);
    winAllSubAllSeedFullCor=tmp.zWinFullCorLasso;
    numROI=156;
    
    indx=load([resultDir, 'zWinAllSeeds_FullCorLasso_',session,'_10min.txt']);
    numClust=length(unique(indx));
    
    disp(['Working on ',session,'_',num2str(numClust),' clusters.'])
    
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
        save([resultDir, 'clusterMean_',num2str(numClust),'clusters_',session,'_normWin.mat'],'finalMeanWinOfClust')
        
        disp('Create brain map of cluster centroids.')
        [Outdata,VoxDim,Header]=rest_readfile([maskDir,'final_reduced.nii']);
        [nDim1 nDim2 nDim3]=size(Outdata);
        temp=unique(Outdata);
        ROIIndx=temp(find(temp~=0));
        numROI=length(ROIIndx);
        
        stateMap=Outdata;
        for m=1:numROI
            stateMap(find(Outdata==ROIIndx(m)))=meanWinOfClust(m);
        end
        
        Header.pinfo = [1;0;0];
        Header.dt    =[16,0];
        rest_WriteNiftiImage(stateMap,Header,[figDir,num2str(numClust),'clusters_state', num2str(i),'_',session,'_normWin.nii']);
    end
end

