clear
clc

session='session1'
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/', session,'/'];
maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/', session,'/thresholdedStatsMap/allSeeds/avgSub/'];

tmp=load([resultDir,'zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'.mat']);
winAllSubAllSeedFullCor=tmp.zWinFullCorLasso;
numROI=156;
numSub=22;
numWinPerSubPerSeed=272;
numSeed=4;
numWinPerSeed=272*4;
numState=5;

indx=load([resultDir, 'clustIndxNormWinAllSeeds_FullCorLasso_',session,'_10min.txt']);
numClust=length(unique(indx));

disp(['Working on ',session,'_',num2str(numClust),' clusters.'])

subAvgiFCWin1=zeros(numSub, numROI);
indxSub=zeros(numWinPerSubPerSeed,numSeed);
winSub=zeros(numWinPerSubPerSeed, numROI, numSeed);
meanWin=zeros(numSub, numROI, numState);
for i=1:numSub
    for j=1:numSeed
        
        indxSub(:,j) =indx(1+numWinPerSubPerSeed*(i-1)+numWinPerSeed*(j-1):numWinPerSubPerSeed*i+numWinPerSeed*(j-1));
        winSub(:, :, j)=winAllSubAllSeedFullCor(1+numWinPerSubPerSeed*(i-1)+numWinPerSeed*(j-1):numWinPerSubPerSeed*i+numWinPerSeed*(j-1), :);
    end
    indxSub1D=reshape(indxSub, [], 1);
    winSub2D=vertcat(winSub(:, :, 1), winSub(:, :, 2), winSub(:, :, 3), winSub(:, :, 4));
    
    for k=1:numState
        indxList=find(indxSub1D==k);
        numWinInClust=length(indxList)
        allWinInClust=zeros(numWinInClust,numROI);
        for n=1:numWinInClust
            for m=1:numROI
                allWinInClust(n,m)=winSub2D(indxList(n), m);
            end
        end
        meanWin(i, :, k)=mean(allWinInClust);
    end
end


for  j=1:numState
    stateMean=zeros(1, numROI);
    pValue=zeros(1, numROI);
    for i=1:numROI
        [h, pTmp]=ttest(meanWin(:, i, j));
        tmp=squeeze(meanWin(:, i, j));
        stateMean(1, i)=nanmean(tmp); % meanDif>0 means state 1 greater than state 2: the order was changed to match the concated data.
        pValue(1, i)=pTmp;
    end
    
    q=0.05;
    [pID,pN] = FDR(pValue,q)
    
    length(find(pValue<=pID))
    
    for k=1:numROI
        p=pValue(k);
        if p<=pID
            if p==0
                p=1e-15;
                negLogPValue(k)=(-1)*log10(p);
                if stateMean(1, k)<0
                    negLogPValue(k) = -negLogPValue(k);
                else
                    negLogPValue(k)=negLogPValue(k);
                end
            else
                negLogPValue(k)=(-1)*log10(p);
                if stateMean(1, k)<0
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
    numSigROI=length(find(negLogPValue~=0))
    
    [Outdata,VoxDim,Header]=rest_readfile([maskDir,'final_reduced.nii']);
    [nDim1 nDim2 nDim3]=size(Outdata);
    temp=unique(Outdata);
    ROIIndx=temp(find(temp~=0));
    numROI=length(ROIIndx);
    
    % plot the -log10(p), the negative numbers were added with a negative
    % sign to display purpose. All -log10(p) were actually positive.
    stateMap=Outdata;
    for m=1:numROI
        stateMap(find(Outdata==ROIIndx(m)))=negLogPValue(m);
    end
    
    Header.pinfo = [1;0;0];
    Header.dt    =[16,0];
    rest_WriteNiftiImage(stateMap,Header,[figDir,'thresholdedStatsMap_state', num2str(j), '_avgSub', session,'_NormWin.nii']);
end





