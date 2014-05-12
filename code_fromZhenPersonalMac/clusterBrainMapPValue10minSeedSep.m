

clear
clc
close all

session='session2';
dataLength='all_10min';
mapType='thresholdedStatsMap';
seedList={'seed1','seed2','seed3','seed4'};
numSeed=length(seedList);
numWinPerSub=272;
numSub=22;
numWinPerSeed=numWinPerSub*numSub;

numROI=156;

for n=1:numSeed
    seed=char(seedList{n});
    disp(['Work on ' seed])
    resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,filesep,session,'/'];
    maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/'];
    figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep,session, filesep, mapType, filesep, seed, '/'];
    
    indx=load([resultDir,'zWinEachSeed_FullCorLasso_',session,'_10min.txt']);
    tmp1=load([resultDir,'zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'.mat']);
    zWinFullCorLasso=tmp1.zWinFullCorLasso;
    
    numClust=length(unique(indx(:,n)));
    disp ('Files loaded successfully.')
    
    for i=1:numClust
        indxList=find(squeeze(indx(:,n))==i);
        numWinInClust=length(indxList);
        allWinInClust=zeros(numWinInClust,numROI);
        for j=1:numWinInClust
            indxClust=indxList(j)+numWinPerSeed*(n-1);
            for m=1:numROI
                allWinInClust(j,m)=zWinFullCorLasso(indxClust, m);
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
        rest_WriteNiftiImage(stateMap,Header,[figDir,'thresholdedStatsMap_',num2str(numClust),'clusters_state', num2str(i),'_',seed, '_', session,'_normWin.nii']);
    end
end



