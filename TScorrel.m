
close all
clear
clc
session={'session2'};
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;

numSub=length(subList);


analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
numSeed=4;

correlAllSub=zeros(numSeed, numSeed, numSub);
for k=1:numSub
    sub=subList{k};
    disp (['Working on sub ', char(sub),' ......'])
    subDir=[analyDir,'data/645/all_10min/', char(session),'/',char(sub)];
    seedROISignals = load([subDir,'/ROISignals_seedROISignal.mat']);
    TC1=seedROISignals.ROISignals;
    correl=corrcoef(TC1);
  correlAllSub(:,:,k)=correl;
end

correlAllSub2D=reshape(correlAllSub,[], numSub)';
meanCorrel=mean(correlAllSub2D);
meanCorrelReshape=reshape(meanCorrel, numSeed, numSeed);
figure(1)
imagesc(meanCorrelReshape)
colorbar
caxis([-0.5 0.5])

a=squeeze(correlAllSub(1,2,:));
b=squeeze(correlAllSub(1,3,:));
c=squeeze(correlAllSub(1,4,:));
d=squeeze(correlAllSub(2,3,:));
e=squeeze(correlAllSub(2,4,:));
f=squeeze(correlAllSub(3,4,:));
correlSeeds=[a,b,c,d,e,f];
save('correlSeedsSession1.txt', '-ascii', '-double', '-tabs','correlSeeds')
    