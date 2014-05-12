function [ numWinInClust ] = countWinInClust(indx)
%This function will compute the number of windows in each cluster
%   Detailed explanation goes here
numClust=length(unique(indx))
numWinInClust=zeros(1,numClust);
for i=1:numClust
    numWinInClust(1,i)=sum(indx==i);
end

