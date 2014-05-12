function [ R2 pseudoF ] = clustNumEstimate10min( session, minNumClust,maxNumClust )
%UNTITLED Summary of this function goes here
%   input: 1. session: string, e.g. session='session1'
% 2. minNumClust: minimum number of clusters, numerical, e.g. minNumClust=2
% 3. maxNumClust: maxmum number of clusters, numerical, e.g.
% maxNumClust=100
% session='session1';
% minNumClust=2;
% maxNumClust=100;
numEstimate=maxNumClust-minNumClust+1;

resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min',filesep, session,'/'];

tmp1=load([resultDir,'zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'.mat']);
zWinFullCorLasso=tmp1.zWinFullCorLasso;

indxList=load([resultDir,'clustIndx_',num2str(minNumClust),'to',num2str(maxNumClust), 'clusters_zWinAllSeeds_FullCorLasso_', session, '_10min.txt']);

disp ('The window and index files are loaded successfully.')

numWin=size(zWinFullCorLasso, 1);
numROI=size(zWinFullCorLasso, 2);

% compute center of all windows
meanzWinFullCorLasso=mean(zWinFullCorLasso);

% compute the sum distance of each window to the ground mean window
T=0;
for i=1:numWin
    for j=1:numROI
        T=T+(zWinFullCorLasso(i,j)-meanzWinFullCorLasso(j))^2;
    end
end

% compute the PG: sum squares within each cluster

PG(1:numEstimate)=0;
R2(1:numEstimate)=0;     %R2 method
pseudoF(1:numEstimate)=0;      %pseudo-F statistic method

for m=minNumClust:maxNumClust
    disp(['Working on ', num2str(m),'clusters.'])
    indx=indxList(:,m-1);
    numClust=m;
    sumWinInClust(1:numClust,1:numROI)=0;  %initial the cluster centers
    numWinInClust(1:numClust)=0;    %initial the number of points in each cluster
    for j=1:numWin
        for k=1:numROI
            sumWinInClust(indx(j,1),k)=sumWinInClust(indx(j,1),k)+zWinFullCorLasso(j,k);   %Summary of all the points which are in the same cluster
        end
        numWinInClust(indx(j))=numWinInClust(indx(j))+1;  %the points number in the cluster
    end
    
    %Wi--the center of the cluster j
    for j=1:numClust
        for k=1:numROI
            sumWinInClust(j,k)=sumWinInClust(j,k)/numWinInClust(j);      %get the center of cluster J in each dimension
        end
    end
    
    %PG--(the SUM distance btween the points to their cluster I's center)^2
    for j=1:numWin
        for k=1:numROI
            PG(m-1)=PG(m-1)+(zWinFullCorLasso(j,k)-sumWinInClust(indx(j,1),k))^2;  %Summary of (distance between points and their cluster center)^2
        end
    end
    
    %R^2 method
    R2(m-1)=1-PG(m-1)/T;
    
    %pseudo-F statistic method
    pseudoF(m-1)=((T-PG(m-1))/(m-1))/(PG(m-1)/(numWin-m));
    
    clear sumWinInclust numWinInClust
end
disp([num2str(numEstimate), ' estimations done!'])

end




