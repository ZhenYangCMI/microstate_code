function [featureWin, zfeatureWin] = featureExtract(FCWinFullBrain)
% This function extract features from the full brain correlation matrix and standardize the featured windows

numWinAllSub=size(FCWinFullBrain, 3);
numROI=size(FCWinFullBrain, 1);
numFeature=numROI*(numROI-1)/2;


featureWin=zeros(numWinAllSub, numFeature);
zfeatureWin=zeros(numWinAllSub, numFeature);
for i=1:numWinAllSub
    disp(['Working on win', num2str(i)])
    
    % full brain correlation of one win
    win=squeeze(FCWinFullBrain(:, :, i));
    
    % extract the lower triangle of the full brain correlation, diagonal
    % was excluded
    n=size(FCWinFullBrain, 1);
    lowerTril=win(find(tril(ones(n,n), -1)));
    featureWin(i,:)=lowerTril;
    % standardize the feature window
    zLowerTril=(lowerTril-repmat(mean(lowerTril), 1, size(lowerTril, 2)))./repmat(std(lowerTril), 1, size(lowerTril, 2));
    
    % save it in the featureWin file
    zfeatureWin(i,:)=zLowerTril;
end
end

