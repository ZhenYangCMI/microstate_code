

clear
clc
close all

session='session2';
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',session,'/'];
maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',session,'/',session,'_ward_euclidean/pBrainMapAvg/'];

tmp2=load([resultDir,'winAllSubAllSeed_partialCor_645_',session,'.mat']);
winAllSubAllSeedPartialCor=tmp2.partWinAllSubAllSeed;
numROI=156;

for k=2:2:16
    numClust=k;
    disp(['Working on ',num2str(numClust),' clusters.'])
    
    if ~exist([resultDir,session,'_ward_euclidean/indx_',num2str(numClust),'clusters.mat'],'file')
        disp ('Running hierarchical clustering')
        distBtwWin= pdist(winAllSubAllSeedPartialCor,'euclidean');
        clustTree = linkage(distBtwWin,'ward');
        indx = cluster(clustTree,'maxclust',numClust,'criterion','distance');
        save([resultDir,'indx_',num2str(numClust),'clusters.mat'],'indx')
    else
        tmp1=load([resultDir,session,'_ward_euclidean/indx_',num2str(numClust),'clusters.mat']);
        indx=tmp1.indx;
        disp ('Index file loaded successfully.')
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
        [h,pValue]=ttest(allWinInClust);
        pValue(find(pValue==0))=1e-300;
        negLogPValue=(-1)*log10(pValue);
        
        %         for m=1:numROI
        %             finalMeanWinOfClust(i,m)=meanWinOfClust(1,m);
        %         end
        %         save([resultDir,'clusterMean_',num2str(numClust),'clusters.mat'],'finalMeanWinOfClust')
        
        [Outdata,VoxDim,Header]=rest_readfile([maskDir,'final_reduced.nii']);
        [nDim1 nDim2 nDim3]=size(Outdata);
        temp=unique(Outdata);
        ROIIndx=temp(find(temp~=0));
        numROI=length(ROIIndx);
        
        stateMap=Outdata;
        for m=1:numROI
            stateMap(find(Outdata==ROIIndx(m)))=negLogPValue(m);
        end
        
        Header.pinfo = [1;0;0];
        Header.dt    =[16,0];
        rest_WriteNiftiImage(stateMap,Header,[figDir,'negLogPValue_',num2str(numClust),'clusters_state', num2str(i),'.nii']);
    end
end


