function [msilhouette] = modified_silhouette(etaMat, solution, solutionVector)
%etaMat is the eta-squared matrix - can be mean, consensus etc. It must be
%the full symmetric matrix (not just lower or upper part)
%solution is a single number corresponding to the current K-level solution, 2-cluster, 3-cluster, 4-cluster etc
%solutionVector is the vector of voxel assigments for the current solution


%written by Clare Kelly September 10th 2009

%%STANDARD SILHOUETTE
%S(i) = (min(AVGD_BETWEEN(i,k)) - AVGD_WITHIN(i)) / max(AVGD_WITHIN(i), min(AVGD_BETWEEN(i,k)))

% where AVGD_WITHIN(i) is the average distance from the i-th point to the
% other points in its own cluster, and AVGD_BETWEEN(i,k) is the average
% distance from the i-th point to points in another cluster k.
% 

%%%MODIFIED SILHOUETTE
TRILMat = tril(etaMat, -1);  %%bottom half of eta-sq matrix
%etaMat is the mean eta-squared matrix across the 38
%subjects

for m=1:solution,  % for each cluster
    eval(['INclusterIndices_' num2str(m) '= find(solutionVector==m);']) %create a variable containing the indices of the voxels in that cluster
end

numSol = [1:solution];

for m=1:solution,

    eval(['currINIndices = INclusterIndices_' num2str(m) ';']); %taking the voxel indices for cluster 1

    tempMat = TRILMat(currINIndices, currINIndices);  %% tempMat contains the eta-sq values for each pairing of voxels
    %within the current cluster, using the tril mat here because do not
    %want i=j, although because I'm taking the mean across all
    %pairings, it doesn't matter
    
    if isnan(sum(sum(tempMat(find(tempMat~=0))))) | (sum(sum(tempMat(find(tempMat~=0)))))==0;,  % all this just takes the mean across all pairings, avoiding any zeros or other problems, and puts
        IN(m) = 0;
    else IN(m) = mean(tempMat(find(tempMat~=0)));
    end

    clear tempMat

    currClust = m;
    otherClust = numSol(find(numSol~=m)); %get all the other clusters

    for oC = 1:length(otherClust), %for each of the other clusters
        eval(['currOUTIndices = INclusterIndices_' num2str(otherClust(oC)) ';']); %get the voxels indices for each other cluster

        tempMat = etaMat(currINIndices, currOUTIndices);  %get eta-sq values for each pairing of voxels. Use full eta-squared matrix because no i = j

        eval(['OUT_cluster(m,' num2str(otherClust(oC)) ') = mean(tempMat(find(tempMat~=0)));']), %get mean of all eta-squared values
        clear tempMat

    end

    msilhouette(m) = (min(1 - OUT_cluster(m, (find(OUT_cluster(m,:)~=0)))) - (1 - IN(m))) / max((1 - IN(m)), min(1 - OUT_cluster(m,(find(OUT_cluster(m,:)~=0)))));
    % This is the same as the standard equation above, but now working
    % on the means for each cluster

    clear OUT*
end

clear IN*

end



