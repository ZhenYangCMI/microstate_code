clear
clc

session='session1'
 resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/', session,'/'];
    maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/'];
    figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/', session,'/DifMapBetwStates/'];
    
    tmp=load([resultDir,'winFullCorLasso_OptimalLambdaPerSub_645_',session,'.mat']);
    winAllSubAllSeedFullCor=tmp.winFullCorLasso;
    numROI=156;
    
    indx=load([resultDir, 'clustIndxNormWinAllSeeds_FullCorLasso_',session,'_10min.txt']);
    numClust=length(unique(indx));
    
    disp(['Working on ',session,'_',num2str(numClust),' clusters.'])
  
    
    % nonParametric test
    for i=1:2
        indxList=find(indx==i);
        numWinInClust=length(indxList);
        allWinInClust=zeros(numWinInClust,numROI);
        for j=1:numWinInClust
            for m=1:numROI
                if i==1
                allWinInClust1(j,m)=winAllSubAllSeedFullCor(indxList(j), m);
                elseif i==2
                  allWinInClust2(j,m)=winAllSubAllSeedFullCor(indxList(j), m);
                end
            end
        end
    end
    
    mean1=mean(allWinInClust1);
    mean2=mean(allWinInClust2);
    Dif=mean2-mean1 % Here if Dif > 0, the results should be interprested as state 1 > state 2, the oder was changed to match the concatednated cluster corder
    
    
       pValue=zeros(1, numROI);     
    for k=1:numROI
        states=mwwtest(allWinInClust1(:,k), allWinInClust2(:,k));
        pValue(1, k)=states.p*2;
    end
    
    q=0.05;
    [pID,pN] = FDR(pValue,q)
    
    length(find(pValue<=pID))
    
        for k=1:numROI
        p=pValue(k);
        if p<=pID
            if p==0
                p=1e-320;
                negLogPValue(k)=(-1)*log10(p);
                if Dif(1,k)<0
                    negLogPValue(k) = -negLogPValue(k);
                else
                    negLogPValue(k)=negLogPValue(k);
                end
            else
                negLogPValue(k)=(-1)*log10(p);
                if Dif(1,k)<0
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
    rest_WriteNiftiImage(stateMap,Header,[figDir,'thresholdedStatsMap_DifState1AndState2_MWWtest', session,'_unNormWin.nii']);
    un
    
    
    
  
  