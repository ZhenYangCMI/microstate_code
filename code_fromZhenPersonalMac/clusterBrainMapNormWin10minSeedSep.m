
clear
clc
close all

sessionList={'session1','session2'};
numSession=length(sessionList)
dataLength='all_10min';
seedList={'seed1','seed2','seed3','seed4'};
numSeed=length(seedList);
numWinPerSub=272;
numSub=22;
numWinPerSeed=numWinPerSub*numSub;


for p=1:numSession
    session=char(sessionList{p});
    disp(['Work on ', session])
    resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength, filesep, session,'/'];
    %     maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/'];
    %     figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep, session,'/corBrainMapAvg/'];
    %
    tmp=load([resultDir,'zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'.mat']);
    winAllSubAllSeedFullCor=tmp.zWinFullCorLasso;
    numROI=156;
    
    indx=load([resultDir, 'zWinEachSeed_FullCorLasso_',session,'_10min.txt']);
    
    for n=1:numSeed
        clear finalMeanWinOfClust
        seed=char(seedList{n});
        disp(['Work on ', seed])
        clustIndxSeed=squeeze(indx(:,n));
        numClust=length(unique(clustIndxSeed))
        
        disp(['Working on ',session,'_',seed, '_',num2str(numClust),' clusters.'])
        
        for i=1:numClust
            indxList=find(clustIndxSeed==i);
            numWinInClust=length(indxList)
            allWinInClust=zeros(numWinInClust,numROI);
            for j=1:numWinInClust
                indxClust=indxList(j)+numWinPerSeed*(n-1);
                for m=1:numROI
                    allWinInClust(j,m)=winAllSubAllSeedFullCor(indxClust, m);
                end
            end
            meanWinOfClust=mean(allWinInClust);
            for m=1:numROI
                finalMeanWinOfClust(i,m)=meanWinOfClust(1,m);
            end
        end
        disp('Cluster centroids computed. Start saving file.')
        save([resultDir, 'clustMean/clusterMean_',num2str(numClust),'clusters_',seed,'_',session,'_normWin.mat'],'finalMeanWinOfClust')
        
        %         disp('Create brain map of cluster centroids.')
        %         [Outdata,VoxDim,Header]=rest_readfile([maskDir,'final_reduced.nii']);
        %         [nDim1 nDim2 nDim3]=size(Outdata);
        %         temp=unique(Outdata);
        %         ROIIndx=temp(find(temp~=0));
        %         numROI=length(ROIIndx);
        %
        %         stateMap=Outdata;
        %         for m=1:numROI
        %             stateMap(find(Outdata==ROIIndx(m)))=meanWinOfClust(m);
        %         end
        %
        %         Header.pinfo = [1;0;0];
        %         Header.dt    =[16,0];
        %         rest_WriteNiftiImage(stateMap,Header,[figDir,num2str(numClust),'clusters_state', num2str(i),'_',session,'_normWin.nii']);
    end
end

