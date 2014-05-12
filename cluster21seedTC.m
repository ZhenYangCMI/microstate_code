clear
clc

numSub=22;
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;
numSub=size(subList,1);
session='session1';
numSeed=21;
numVol=884;

covType='GSR'

analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/', covType, filesep];
resultsDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/', covType, filesep];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/', covType, '/', session, '/21seedTS/'];
TCAllSubCat=zeros(numVol,numSeed);
TCAllSub=zeros(numVol, numSeed, numSub);
for i=1:numSub
    sub=subList{i};
    disp (['Working on sub ', char(sub),' ......'])
    subDir=[analyDir, session,'/',char(sub)];
    seedROISignals = load([subDir,'/ROISignals_21seedROISignal.mat']);
    TC=seedROISignals.ROISignals;
    numSeed=size(TC,2);
    TCAllSubCat=vertcat(TCAllSubCat,TC);
    TCAllSub(:,:,i)=TC;
end

TCAllSub2D=reshape(TCAllSub, [],numSub)';
meanTCAllSub=mean(TCAllSub2D);
mean2D=reshape(meanTCAllSub, numVol, numSeed);

% clust the normalized TS
normMean2D=(mean2D-repmat(mean(mean2D), size(mean2D, 1), 1))./repmat(std(mean2D), size(mean2D, 1), 1);
normMean2D=normMean2D';

%save([resultsDir, session, '/21seedTS/normMeanTS21seed.mat'], 'normMean2D')

disp('run hierachical clustering')
D = pdist(normMean2D,'euclidean');
disp('Distance computed')
Z=linkage(D, 'ward');
figure(1)
dendrogram(Z)
T = cluster(Z,'maxclust',5)
disp('cluster index generated')
saveas(figure(1), [figDir, 'zMeanTS21Seed_', session,'_', covType, '_matlab.png'])

% cluster the unnormalized TS
% mean2D=mean2D';
% disp('run hierachical clustering')
% D = pdist(mean2D,'euclidean');
% disp('Distance computed')
% Z=linkage(D, 'ward');
% dendrogram(Z)
% T = cluster(Z,'maxclust',4)
% disp('cluster index generated')


