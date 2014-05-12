clear
clc

% load subList
subList={'3808535','8574662'};
numSub=length(subList);
numROI='0100';

session='session1'
winType='winFullCorLasso'

% define the windodw parameters
% below for TR645
% winSize=69;
% step=3; % in TR
% numVol=884;

winSize=18;
step=1; % in TR
numVol=116;

analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/'];
dataDir=[analyDir, 'testing/TR645_5min/'];
resultDir=[analyDir,'testing/results/TR645_5min/'];

tmp=load([resultDir, sprintf('featureFC_%s_%s_%s.mat', winType, numROI, session)]) % unstandarded feature windows were used following Allen et al.
winFullCorLasso=tmp.featureWin;
variance=var(winFullCorLasso');
[pks locs]=findpeaks(variance);

figure(1)
plot(variance)

examplar=zeros(length(locs), size(winFullCorLasso,2));
for i=1:length(locs)
    winIndx=locs(i);
    examplar(i,:)=winFullCorLasso(winIndx, :);
end
minNumClust=1;
maxNumClust=20;
%CVI=zeros(maxNumClust, 1);
index=zeros(size(winFullCorLasso,1), (maxNumClust-minNumClust+1));

for numClust=minNumClust:maxNumClust
    clear indx ctrs_stat sumd D
    opts=statset('MaxIter',500);
    numReplications=500
    [indx1,ctrs1]=kmeans(examplar,numClust,'Distance','cityblock','emptyaction','singleton','Replicates',numReplications,'options',opts);
    initialCtrs=repmat(ctrs1, [1,1,numReplications]);
    [indx2,ctrs2,sumd,D]=kmeans(winFullCorLasso,numClust,'Distance','cityblock','emptyaction','singleton','Replicates',numReplications,'options',opts, 'start', initialCtrs);
    % compute the cluster validity index
    CVI(numClust,1)=sum(sumd)./(sum(sum(D))-sum(sumd));
    index(:,numClust)=indx2;
    disp (['k-means clustering with', numROI, '_', num2str(numClust),' cluster done.'])
end
logCVI=log10(CVI);
save([resultDir, sprintf('CVI_%ito%iclusters_%ireplications_%sROIs.mat', minNumClust, maxNumClust, numReplications, numROI)], 'CVI')
save([resultDir, sprintf('logCVI_%ito%iclusters_%ireplications_%sROIs.mat', minNumClust, maxNumClust, numReplications, numROI)], 'logCVI')
save([resultDir, sprintf('clustIndx_%ito%iclusters_%ireplications_%sROIs.mat', minNumClust, maxNumClust, numReplications, numROI)], 'index')

figure(2)
%subplot(1,2,1)
plot(CVI)
% subplot(1,2,2)
% plot(logCVI)
saveas(figure(2), [resultDir, sprintf('CVI_%ito%iclusters_%ireplications_%sROIs.png', minNumClust, maxNumClust, numReplications, numROI)])


