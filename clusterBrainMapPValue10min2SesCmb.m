

clear
clc
close all

session='2sessions';
dataLength='all_10min';
dataType='even'
mapType='thresholdedStatsMap';
seedNum='allSeeds';

numROI=156;
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,filesep,session,'/'];
maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep,session, filesep, mapType, filesep];

indx=load([resultDir,'zWinAllSeeds_FullCorLasso_',session,'_10min_', dataType, '.txt']);
tmp1=load([resultDir,'zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'_', dataType,'.mat']);
zWinFullCorLasso=tmp1.(dataType);
numClust=length(unique(indx))
disp ('Files loaded successfully.')

for i=1:numClust
    indxList=find(indx==i);
    numWinInClust=length(indxList);
    allWinInClust=zeros(numWinInClust,numROI);
    for j=1:numWinInClust
        for m=1:numROI
            allWinInClust(j,m)=zWinFullCorLasso(indxList(j), m);
        end
    end
    meanWinOfClust=mean(allWinInClust);
    [h,pValue]=ttest(allWinInClust);
    [pID,pN] = FDR(pValue,0.05)
    
    for k=1:numROI
        p=pValue(k);
        if p<=pID
            if p==0
                p=1e-320;
                negLogPValue(k)=(-1)*log10(p);
                if meanWinOfClust(k)<0
                    negLogPValue(k) = -negLogPValue(k);
                else
                    negLogPValue(k)=negLogPValue(k);
                end
            else
                negLogPValue(k)=(-1)*log10(p);
                if meanWinOfClust(k)<0
                    negLogPValue(k) = -negLogPValue(k);
                else
                    negLogPValue(k)=negLogPValue(k);
                end
            end
        else
            negLogPValue(k)=0;
        end
    end
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
    rest_WriteNiftiImage(stateMap,Header,[figDir,'thresholdedStatsMap_',num2str(numClust),...
        'clusters_state', num2str(i),'_',session,'_', dataType, '_normWin.nii']);
end



