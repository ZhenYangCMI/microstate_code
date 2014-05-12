
clear
clc
close all

session='session1';
lambda=0.12;
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',session,'/'];
maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',session,'/',session,'_ward_euclidean/'];

tmp2=load([resultDir,'winAllSubAllSeed_partialCor_645_',session,'_',num2str(lambda),'.mat']);
winAllSubAllSeedPartialCor=tmp2.winAllSubAllSeedPartialCor;
numROI=156;

minNumClust=2;
maxNumClust=2;
for k=minNumClust:maxNumClust
    numClust=k;
    disp(['Working on ',session,'_',num2str(numClust),' clusters.'])
    if ~exist([resultDir,session,'_ward_euclidean/lambda',num2str(lambda),'/indxList_',num2str(minNumClust),'_to_',num2str(maxNumClust),'_clusters_',num2str(lambda),'.mat'],'file')
        disp('Compute indxList')
        [ indxList clustTreeList ] = computeClustIndx( minNumClust, maxNumClust, session );
    else
        tmp1=load([resultDir,session,'_ward_euclidean/lambda',num2str(lambda),'/indxList_',num2str(minNumClust),'_to_',num2str(maxNumClust),'_clusters_',num2str(lambda),'.mat']);
        indx=tmp1.indxList;
        disp('indxList loaded successfully.')
    end
    
    for i=1:numClust
        indxList=find(indx==i);
        numWinInClust=length(indxList);
        allWinInClust=zeros(numWinInClust,numROI);
        for j=1:numWinInClust
            for m=1:numROI
                allWinInClust(j,m)=winAllSubAllSeedPartialCor(indxList(j), m);
            end
        end
        meanWinOfClust=mean(allWinInClust);
        for m=1:numROI
            finalMeanWinOfClust(i,m)=meanWinOfClust(1,m);
        end
        disp('Cluster centroids computed. Start saving file.')
        save([resultDir,session,'_ward_euclidean/lambda',num2str(lambda),'/clusterMean_',num2str(numClust),'clusters.mat'],'finalMeanWinOfClust')
        
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
        rest_WriteNiftiImage(stateMap,Header,[figDir,num2str(numClust),'clusters_state', num2str(i),'.nii']);
    end
end

