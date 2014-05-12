function [ indxList clustTreeList ] = computeClustIndx( minNumClust, maxNumClust, session )
% this function output an indx maxtrix (numWin * maxNumClust)and
% clustTree (numWin-1)*3 * maxNumclust
% input: 1. minNumClust: numerical, e.g minNumClust=1
%2. maxNumClust: numerical, e.g. maxNumClust=100;
% 3. session: string, e.g. session='session1';
% minNumclust=1
% maxNumClust=5
% session='session1'

if strcmp(session,'session1')
    lambda=0.12
else
    lambda=0.11
end
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',session,'/'];

tmp1=load([resultDir,'winAllSubAllSeed_partialCor_645_',session,'_',num2str(lambda),'.mat']);
winAllSubAllSeedPartialCor=tmp1.winAllSubAllSeedPartialCor;
numWin=size(winAllSubAllSeedPartialCor, 1);
numIter=maxNumClust-minNumClust+1;
clustTreeList=zeros(numWin-1, 3, numIter);
indxList=zeros(numWin,numIter);

k=0;
for i=minNumClust:maxNumClust
    numClust=i;
    k=k+1
    disp(['Working on ',num2str(numClust), ' clusters.'])
    distBtwWin= pdist(winAllSubAllSeedPartialCor,'euclidean');
    clustTree = linkage(distBtwWin,'ward');
    indx = cluster(clustTree,'maxclust',numClust,'criterion','distance');
    for j=1:numWin-1
        for m=1:3
            clustTreeList(j,m,k)=clustTree(j,m);
        end
    end
    for j=1:numWin
        indxList(j,k)=indx(j,1);
    end
end
disp ('Clustering done. Start saving results ...')
save([resultDir,session,'_ward_euclidean/lambda',num2str(lambda),'/','indxList_',num2str(minNumClust),'_to_', num2str(maxNumClust),'_clusters_',num2str(lambda),'.mat'],'indxList')
save([resultDir,session,'_ward_euclidean/lambda',num2str(lambda),'/','clustTreeList_',num2str(minNumClust),'_to_', num2str(maxNumClust),'_clusters_',num2str(lambda),'.mat'],'clustTreeList')
disp (['indx and clustTree saved under ', resultDir, session, '_ward_euclidean'])