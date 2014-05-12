

clear
clc
close all

session='session2';
dataLength='all_10min';
%mapType='thresholdedStatsMap';
mapType='thresholdedStatsMapR'
seedNum='allSeeds';

numROI=156;
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,filesep,session,'/'];
maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep,session, filesep, mapType, filesep];

indx=load([resultDir,'clustIndxNormWinAllSeeds_FullCorLasso_',session,'_10min.txt']);
tmp1=load([resultDir,'zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'.mat']);
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
    
    % plot the -log10(p), the negative numbers were added with a negative
    % sign to display purpose. All -log10(p) were actually positive.
%     stateMap=Outdata;
%     for m=1:numROI
%         stateMap(find(Outdata==ROIIndx(m)))=negLogPValue(m);
%     end
%     
%     Header.pinfo = [1;0;0];
%     Header.dt    =[16,0];
%     rest_WriteNiftiImage(stateMap,Header,[figDir,'thresholdedStatsMap_',num2str(numClust),'clusters_state', num2str(i),'_',seedNum, '_', session,'_normWin.nii']);
%     
    
    % plot the rank of -log10(p)
    negLogPValueR=zeros(1, numROI);
     negLogPValueR(find(negLogPValue>0))=tiedrank(negLogPValue(find(negLogPValue>0)));
     negLogPValueR(find(negLogPValue<0))=-tiedrank(-negLogPValue(find(negLogPValue<0)));
    stateMapR=Outdata;
    for m=1:numROI
        stateMapR(find(Outdata==ROIIndx(m)))=negLogPValueR(m);
    end
    max(negLogPValueR)
    Header.pinfo = [1;0;0];
    Header.dt    =[16,0];
    rest_WriteNiftiImage(stateMapR,Header,[figDir,'thresholdedStatsMapR_',num2str(numClust),'clusters_state', num2str(i),'_',seedNum, '_', session,'_normWin.nii']);
        
end



