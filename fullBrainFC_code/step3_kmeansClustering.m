clear
clc

% load subList
clear
clc
close all

% load subList
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;
%subList={'0021002', '0021006'};
numSub=length(subList);

session='session1';

% define the analysis strategy, use to fid the right TS data
covType='GSR';


% define the windodw parameters
% below for TR645
winSize=69;
step=3; % in TR
numVol=884;


analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/'];
resultDir=[analyDir,'results/', covType, '/', session, '/'];

tmp=load([resultDir, sprintf('featureFC_winFullCorLasso_%s.mat',session)])
winFullCorLasso=tmp.featureWin;
variance=var(winFullCorLasso');
[pks locs]=findpeaks(variance);

figure(1)
plot(variance)

% extract the examplar windows
examplar=zeros(length(locs), size(winFullCorLasso,2));
for i=1:length(locs)
    winIndx=locs(i);
    examplar(i,:)=winFullCorLasso(winIndx, :);
end

% kmeans clustering the feature windows in two steps
minNumClust=1;
maxNumClust=20;
CVI=zeros(maxNumClust, 1);
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
       disp (['k-means clustering with', num2str(numClust),' cluster done.'])
end
logCVI=log10(CVI);
save([resultDir, sprintf('CVI_%ito%iclusters_%ireplications_%s.mat', minNumClust, maxNumClust, numReplications, session)], 'CVI')
save([resultDir, sprintf('logCVI_%ito%iclusters_%ireplications_%s.mat', minNumClust, maxNumClust, numReplications, session)], 'logCVI')
save([resultDir, sprintf('clustIndx_%ito%iclusters_%ireplications_%s.mat', minNumClust, maxNumClust, numReplications, session)], 'index')

figure(2)
subplot(1,2,1)
plot(CVI)
subplot(1,2,2)
plot(logCVI)
saveas(figure(2), sprintf('CVI_%ito%iclusters_%ireplications_%s.png', minNumClust, maxNumClust, numReplications, session))


