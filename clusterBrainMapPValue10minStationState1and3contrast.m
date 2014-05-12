

clear
clc
close all

session='session1';
dataLength='all_10min';
%mapType='stationaryFCMap';
mapType='stationaryFCMap';

numSeed=4;
numROI=156;
numSub=22;
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,filesep,session,'/'];
maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep,session, filesep, mapType, filesep];

tmp=load([resultDir,'/z_stationary_FC.mat']);
stationaryFC=tmp.z_stationary_FC;
disp ('Files loaded successfully.')


stationaryFCSeed1=zeros(numSub,numROI);
stationaryFCSeed3=zeros(numSub,numROI);
for j=1:numSub
    for m=1:numROI
        stationaryFCSeed1(j,m)=stationaryFC(1,m,j);
        stationaryFCSeed3(j,m)=stationaryFC(3,m,j);
    end
end
difFC=stationaryFCSeed1-stationaryFCSeed3;
[h,pValue, ci, stats]=ttest(stationaryFCSeed1,stationaryFCSeed3);
[pID,pN] = FDR(pValue,0.05)

t=stats.tstat;
for k=1:numROI
    p=pValue(k);
    if p<=pID
       t(k)=t(k);
                    else
        t(k)=0;
    end
end
max(t)
min(t)

[Outdata,VoxDim,Header]=rest_readfile([maskDir,'final_reduced.nii']);
[nDim1 nDim2 nDim3]=size(Outdata);
temp=unique(Outdata);
ROIIndx=temp(find(temp~=0));
numROI=length(ROIIndx);
%
%     % plot the rank of -log10(p)
%     negLogPValueR=zeros(1, numROI);
%      negLogPValueR(find(negLogPValue>0))=tiedrank(negLogPValue(find(negLogPValue>0)));
%      negLogPValueR(find(negLogPValue<0))=-tiedrank(-negLogPValue(find(negLogPValue<0)));
%
%      numPos=length((find(negLogPValue>0)))
%      numNeg=length((find(negLogPValue>0)))
% plot rank of the -log10(p) map
%      stationaryFCMapR=Outdata;
%     for m=1:numROI
%         stationaryFCMapR(find(Outdata==ROIIndx(m)))=negLogPValueR(m);
%     end
%     max(abs(negLogPValueR))
%     Header.pinfo = [1;0;0];
%     Header.dt    =[16,0];
%rest_WriteNiftiImage(stationaryFCMapR,Header,[figDir,'thresholdedStationaryFCMapR_seed',num2str(i),'_', session, '.nii']);

% plot -log10(p) map
stationaryFCMap=Outdata;
for m=1:numROI
    stationaryFCMap(find(Outdata==ROIIndx(m)))=t(m);
end

Header.pinfo = [1;0;0];
Header.dt    =[16,0];
%rest_WriteNiftiImage(stationaryFCMap,Header,[figDir,'thresholdedStationaryFCMap_seed1-seed3_', session, '.nii']);

rest_WriteNiftiImage(stationaryFCMap,Header,['thresholdedStationaryFCMap_seed1-seed3_', session, '.nii']);



