
clear
clc
close all

sessionList={'2sessions'};
numSession=length(sessionList)
dataLength='all_10min';
%dataTypeList={'odd','even'};
dataTypeList={'21sub'};
varName='zWinFullCorLasso21Sub';
winSize=136;
numDataType=length(dataTypeList);
for i=1:numSession
    session=char(sessionList{i});
    resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength, filesep, session,'/'];
    maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/'];
    figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep, session,'/corBrainMapAvg/'];
    for j=1:numDataType
        dataType=char(dataTypeList{j});
        tmp=load([resultDir,'zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'_',dataType,'_win', num2str(winSize),'.mat']);
        %winAllSubAllSeedFullCor=tmp.(dataType);
                winAllSubAllSeedFullCor=tmp.(varName);
        numROI=156;
        
        indx=load([resultDir, 'clustIndxNormWinAllSeeds_FullCorLasso_',session,'_10min_',dataType,'_win', num2str(winSize), '.txt']);
        numClust=length(unique(indx));
        
        disp(['Working on ',session,'_',dataType, '_',num2str(numClust),' clusters.'])
        
        for k=1:numClust
            indxList=find(indx==k);
            numWinInClust=length(indxList);
            allWinInClust=zeros(numWinInClust,numROI);
            for n=1:numWinInClust
                
                allWinInClust(n,:)=winAllSubAllSeedFullCor(indxList(n), :);
                
            end
            meanWinOfClust=mean(allWinInClust);
            for m=1:numROI
                finalMeanWinOfClust(k,m)=meanWinOfClust(1,m);
            end
            
            disp('Create brain map of cluster centroids.')
            %             [Outdata,VoxDim,Header]=rest_readfile([maskDir,'final_reduced.nii']);
            %             [nDim1 nDim2 nDim3]=size(Outdata);
            %             temp=unique(Outdata);
            %             ROIIndx=temp(find(temp~=0));
            %             numROI=length(ROIIndx);
            %
            %             stateMap=Outdata;
            %             for m=1:numROI
            %                 stateMap(find(Outdata==ROIIndx(m)))=meanWinOfClust(m);
            %             end
            %
            %             Header.pinfo = [1;0;0];
            %             Header.dt    =[16,0];
            %             rest_WriteNiftiImage(stateMap,Header,[figDir,num2str(numClust),'clusters_state', num2str(k),'_',session,'_', dataType, '_normWin.nii']);
            %         end
            %         disp('Cluster centroids computed. Start saving file.')
            save([resultDir, sprintf('clusterMean_%dclusters_%s_%s_win%d_normWin.mat',numClust,session,dataType,winSize)],'finalMeanWinOfClust')
        end
    end
end
