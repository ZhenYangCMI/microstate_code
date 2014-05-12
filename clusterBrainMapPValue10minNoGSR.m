

clear
clc
close all

session='2sessions';
dataLength='all_10min';
seedNum='allSeeds';
covType='noGSR'; %'noGSR', 'compCor'

winSize=69;

numROI=156;
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,'/', covType, '/', session,'/'];
maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, '/', covType, '/', session, '/'];

indx=load([resultDir,'clustIndxNormWinAllSeeds_FullCorLasso_',session,'win', num2str(winSize),'_', covType, '.txt']);
tmp1=load([resultDir,'zwinFullCorLasso_OptimalLambdaPerSub_645_',session,'win',num2str(winSize),'_', covType, '.mat']);
zWinFullCorLasso=tmp1.zWinFullCorLasso;
numClust=length(unique(indx));
disp ('Files loaded successfully.')

finalMeanWinOfClust=zeros(numClust, numROI);
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
    finalMeanWinOfClust(i, :)=meanWinOfClust;
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
    rest_WriteNiftiImage(stateMap,Header,[figDir,sprintf('thresholdedStatsMap_%dclusters_state%d_%s_%s_normWin_win%d.nii',numClust,i,covType,session,winSize)]);
        
end
save([resultDir, sprintf('clusterMean_%dclusters_%s_%s_win%d_normWin.mat',numClust,session,covType,winSize)],'finalMeanWinOfClust')



