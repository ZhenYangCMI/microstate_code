

clear
clc
close all

session='session1';
dataLength='all_10min';
mapType='thresholdedStatsMap';
seedNum='allSeeds';

winSize=136;

numROI=156;
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,filesep,session,'/'];
maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep,session, filesep, mapType, filesep];

indx=load([resultDir,'clustIndxWinAllSeeds_FullCorLasso_',session,'_10min_win', num2str(winSize),'.txt']);
tmp1=load([resultDir,'zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'win',num2str(winSize),'.mat']);
zWinFullCorLasso=tmp1.zWinFullCorLasso;
numClust=length(unique(indx));
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
    pValue;
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
    negLogPValue;
    numSigROI1=length(find(negLogPValue~=0))
    
    [Outdata,VoxDim,Header]=rest_readfile([maskDir,'final_reduced.nii']);
    [nDim1 nDim2 nDim3]=size(Outdata);
    temp=unique(Outdata);
    ROIIndx=temp(find(temp~=0));
    numROI=length(ROIIndx);
    
%     plot the -log10(p), the negative numbers were added with a negative
%     sign to display purpose. All -log10(p) were actually positive.
    stateMap=Outdata;
    for m=1:numROI
        stateMap(find(Outdata==ROIIndx(m)))=negLogPValue(m);
    end
    
    Header.pinfo = [1;0;0];
    Header.dt    =[16,0];
    rest_WriteNiftiImage(stateMap,Header,[figDir,sprintf('thresholdedStatsMap_%dclusters_state%d_%s_%s_normWin_win%d.nii',numClust,i,seedNum,session,winSize)]);
        
end



