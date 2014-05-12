%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script do the following
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

% load subList
subList={'3808535','8574662'};
numSub=length(subList);
ROI='0200';
session='session1'
TR='645';
mask='/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/mask/k0200.nii';

close all
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/testing/results/TR', TR, '_5min/']
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/testing/figs/'];
% load the index file

clustIndx1=load([resultDir, 'clustIndx_winFullCorLasso_', ROI, '_2to14clusters.txt']);

clustIndx2=load([resultDir, 'clustIndx_1to20clusters_500replications_', ROI, 'ROIs.mat']);
clustIndx2=clustIndx2.index;

win=load([resultDir, 'zfeatureFC_winFullCorLasso_', ROI, '_', session, '.mat'])
winData=win.zfeatureWin;
numFeature=size(winData, 2);

minNumClust=7
maxNumClust=14
n=maxNumClust-minNumClust+1;
t=0;

for i=minNumClust:maxNumClust
    disp('hierachical')
    numClust=i
       
    modul1=zeros(n, 189);
        for j=1:numClust
        clust1=winData(find(clustIndx1(:,i)==j), :);
        clustMedian1=squareform(median(clust1));
        [Ci Q] = modularity_louvain_und_sign(clustMedian1);
        
        [Outdata,VoxDim,Header]=rest_readfile(mask);
        [nDim1 nDim2 nDim3]=size(Outdata);
        temp=unique(Outdata);
        ROIIndx=temp(find(temp~=0));
        numROI=length(ROIIndx);
        
        modul1Map=Outdata;
        for m=1:numROI
            modul1Map(find(Outdata==ROIIndx(m)))=Ci(m);
        end
        
        Header.pinfo = [1;0;0];
        Header.dt    =[16,0];
        rest_WriteNiftiImage(modul1Map,Header,[figDir,sprintf('hierarchical_statModulMap_%dclusters_state%d_%sTR_%sROIs.nii',i,j,TR, ROI)]);
        
    end
end

for i=minNumClust:maxNumClust
    disp('kmeans')
    numClust=i
       
    modul2=zeros(n, 189);
        for j=1:numClust
        clust2=winData(find(clustIndx2(:,i)==j), :);
        clustMedian2=squareform(median(clust2));
        [Ci Q] = modularity_louvain_und_sign(clustMedian2);
        
        [Outdata,VoxDim,Header]=rest_readfile(mask);
        [nDim1 nDim2 nDim3]=size(Outdata);
        temp=unique(Outdata);
        ROIIndx=temp(find(temp~=0));
        numROI=length(ROIIndx);
        
        modul2Map=Outdata;
        for m=1:numROI
            modul2Map(find(Outdata==ROIIndx(m)))=Ci(m);
        end
        
        Header.pinfo = [1;0;0];
        Header.dt    =[16,0];
        rest_WriteNiftiImage(modul2Map,Header,[figDir,sprintf('kmeans_statModulMap_%dclusters_state%d_%sTR_%sROIs.nii',i,j,TR, ROI)]);
        
    end
end




