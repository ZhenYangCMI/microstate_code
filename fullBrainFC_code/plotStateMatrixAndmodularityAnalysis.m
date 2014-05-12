%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script do the following
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

% load subList
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;
resultDir='/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/results/noGSRFiltered/';
mask='/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/mask/AAL_61x73x61_YCG_reduced.nii'
figDir = '/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/figs/noGSRFiltered/';

clustIndx=load([resultDir, 'clustIndx_normFeature2_session1.txt']);
win=load([resultDir, 'normFeature2_session1.mat'])
winData=win.normFeature2;
numFeature=size(winData, 2);
numClust=length(unique(clustIndx))

for i=1:numClust
    clust=winData(find(clustIndx==i), :);
    clustMedian=squareform(median(clust));
    
%     figure(i)
%     imagesc(clustMedian)
%     colorbar
%     saveas(figure(i), [figDir, 'clustMedian_clust', num2str(i), '.jpg'])

%     if i==1
%     caxis([-0.4 0.4])
%     elseif i==2
%         caxis([-0.4 0.4])
%     elseif i==3
%         caxis([-0.6 0.2])
%     else
%         caxis([0 2])
%     end
        
    [Ci Q] = modularity_louvain_und_sign(clustMedian);
    max(Ci)
%     
%     [Outdata,VoxDim,Header]=rest_readfile(mask);
%     [nDim1 nDim2 nDim3]=size(Outdata);
%     temp=unique(Outdata);
%     ROIIndx=temp(find(temp~=0));
%     numROI=length(ROIIndx);
%     
%     modulMap=Outdata;
%     for m=1:numROI
%         modulMap(find(Outdata==ROIIndx(m)))=Ci(m);
%     end
%     
%     Header.pinfo = [1;0;0];
%     Header.dt    =[16,0];
%     rest_WriteNiftiImage(modulMap,Header,[figDir,sprintf('hierarchical_statModulMap_%dclusters_state%d.nii',numClust, i)]);
    
end



