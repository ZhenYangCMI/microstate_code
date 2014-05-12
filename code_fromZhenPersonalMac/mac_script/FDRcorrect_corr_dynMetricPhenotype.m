clear
clc


numSeed=4;
numPair=6;
numSub=21;
numSubExclude=5;
numState=5;
% sub 1-4 and 11 were removed, after removing each one the order will be different, so the Indx is different from the orignal order
subIDList=[1,1,1,1,7];
resultDir='/Users/zhenyang/Desktop/Zhen/results/all_10min/2sesConcate/';
NPList={'DKEFS_DESIGN_SWITCHING_SCALED', 'DKEFS_VERB_CATEGORYSWITCH_SCALED', 'DKEFS_TOTALWEIGHTACHIEVESCORE_SCALED',...
    'DKEFS_COLORWORD_INHIBITION_SCALED','WASI_PERF','WASI_VERB','WASI_FULL_4'};
dynMetricList1={'transitions', 'TCPrctOverlapStateSum'};
dynMetricList2={'subSeedPrct', 'stateFreq', 'duration'};
numNP=length(NPList);
num2DMetric=length(dynMetricList1);
num3DMetric=length(dynMetricList2);
NPdata1=load([resultDir,'DKEFS_WASI_16sub_2sesCmb_session1.txt']);
NPdata2=load([resultDir,'DKEFS_WASI_16sub_2sesCmb_session2.txt']);


% compute Pval on metric with 2D structure
ICCREML2D=zeros(1, numSeed+numPair);
t=0;
for j=1:num2DMetric
    measure=char(dynMetricList1{j});
    if strcmp(measure,'TCPrctOverlapStateSum')
        numLoop=numPair; %numSeed=numSeedPairs
    else
        numLoop=numSeed;
    end
    % compute ICC and perform ttest on % time spent in a state
    tmp1=load([resultDir, measure,'_session1.mat']);
    data1=tmp1.(measure);
    tmp2=load([resultDir, measure,'_session2.mat']);
    data2=tmp2.(measure);
    data16sub1=data1;
    data16sub2=data2;
    for p=1:numSubExclude
        subID=subIDList(p);
        data16sub1(subID,:)=[];
        data16sub2(subID,:)=[];
    end
    
        for i=1:numLoop
        disp(['Work on seed ', num2str(i)])
        session1Data=squeeze(data1(:,i));
        session2Data=squeeze(data2(:,i));
        Y=[session1Data;session2Data];
        time = [ones(numSub,1);2*ones(numSub,1)];
        sID=[[1:numSub]';[1:numSub]'];
        [ ICC, idx_fail] = do_ICC(Y, time, [], [], sID);
        ICCREML2D(i,j)=ICC;
        
            for m=1:numNP
                
                [RHO1,PVAL1] = corr(NPdata1(:,m),data16sub1(:,i),'Type','Pearson', 'rows','pairwise');
                [RHO2,PVAL2] = corr(NPdata2(:,m),data16sub2(:,i),'Type','Pearson','rows','pairwise');
                PVAL2Dses1(m,i,j)=PVAL1;
                PVAL2Dses2(m,i,j)=PVAL2;
            end
        
    end
end
PVAL2Dses1
PVAL2Dses2

for k=1:num2DMetric
    for i=1:numPair
                    if ICCREML2D(i,j)<0.4
                PVAL2Dses1(:,i,j)=0;
PVAL2Dses2(:,i,j)=0;
            end
        end
end

PVAL2Dses1
PVAL2Dses2
numZerosSes1=length(find(PVAL2Dses1==0))
numZerosSes2=length(find(PVAL2Dses2==0))
numP1=num2DMetric*numPair*numNP-numZerosSes1
numP2=num2DMetric*numPair*numNP-numZerosSes2

PVAL2Dses1(find(PVAL2Dses1==0))=[];
PVAL2Dses2(find(PVAL2Dses2==0))=[];
numPVALfinal=length(PVAL2Dses1)
numPVALfinal=length(PVAL2Dses2)    


% compute Pval on metric with 3D structure

ICCREML3D=zeros(numState,numSeed, num3DMetric);
for k=1:num3DMetric
    measure=char(dynMetricList2{k});
    
    tmp1=load([resultDir, measure,'_session1.mat']);
    data1=tmp1.(measure);
    tmp2=load([resultDir, measure,'_session2.mat']);
    data2=tmp2.(measure);
    data16sub1=data1;
    data16sub2=data2;
    
    for p=1:numSubExclude
        subID=subIDList(p);
        data16sub1(:,:,subID)=[];
        data16sub2(:,:,subID,:)=[];
    end
    
    for i=1:numSeed
        disp(['Work on seed ', num2str(i)])
        for j=1:numState
            session1Data=squeeze(data1(j,i,:));
            session2Data=squeeze(data2(j,i,:));
            Y=[session1Data;session2Data];
            time = [ones(numSub,1);2*ones(numSub,1)];
            sID=[[1:numSub]';[1:numSub]'];
            [ ICC1, idx_fail] = do_ICC(Y, time, [], [], sID);
            ICCREML3D(j,i,k)=ICC1;
            for m=1:numNP
                [RHO1,PVAL3] = corr(NPdata1(:,m),squeeze(data16sub1(j,i,:)),'Type','Pearson', 'rows','pairwise');
                [RHO2,PVAL4] = corr(NPdata2(:,m),squeeze(data16sub2(j,i,:)),'Type','Pearson','rows','pairwise');
                PVAL3Dses1(m,j,i,k)=PVAL3;
                PVAL3Dses2(m,j,i,k)=PVAL4;
            end
        end
    end
end
PVAL3Dses1
PVAL3Dses2


for k=1:num3DMetric
    for i=1:numSeed
        for j=1:numState
            if ICCREML3D(j,i,k)<0.4
                PVAL3Dses1(:,j,i,k)=0;
PVAL3Dses2(:,j,i,k)=0;
            end
        end
    end
end

PVAL3Dses1
PVAL3Dses2
numZerosSes1=length(find(PVAL3Dses1==0))
numZerosSes2=length(find(PVAL3Dses2==0))
numP1=num3DMetric*numSeed*numState*numNP-numZerosSes1
numP2=num3DMetric*numSeed*numState*numNP-numZerosSes2

PVAL3Dses1(find(PVAL3Dses1==0))=[];
PVAL3Dses2(find(PVAL3Dses2==0))=[];
numPVALfinal=length(PVAL3Dses1)
numPVALfinal=length(PVAL3Dses2)

% concate PVALUES
PVALFinal1=[PVAL2Dses1, PVAL3Dses1];
PVALFinal2=[PVAL2Dses2, PVAL3Dses2];

[pID1,pN1] = FDR(PVALFinal1,0.05)
[pID2,pN2] = FDR(PVALFinal2,0.05)

find(PVALFinal1<=pID1)
find(PVALFinal2<=pID2)
PVALses1Survive=PVALFinal1(find(PVALFinal1<=pID1))
PVALses2Survive=PVALFinal2(find(PVALFinal2<=pID2))