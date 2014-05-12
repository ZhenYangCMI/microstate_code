

clear
clc
close all

session='session1';
dataLength='first_5min';
numROI=156;
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,filesep,session,'/'];
maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength,filesep, session,'/',session,'_ward_euclidean/lambdaOptimalPerSub/thresholdedStatsMap/'];

indx=load([resultDir,session,'_ward_euclidean/lambdaOptimalPerSub/winAllSeeds_FullCorLasso_',session,'.txt']);
tmp1=load([resultDir,'winFullCorLasso_OptimalLambdaPerSub_645_',session,'.mat']);
zWinFullCorLasso=tmp1.winFullCorLasso;
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
            else
                negLogPValue(k)=(-1)*log10(p);
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
        rest_WriteNiftiImage(stateMap,Header,[figDir,'thresholdedStatsMap_',num2str(numClust),'clusters_state', num2str(i),'_',session,'.nii']);
    end



