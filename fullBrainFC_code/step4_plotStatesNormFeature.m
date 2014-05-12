%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script do the following for each session separately:
%1. Extract three types of Fisher Z transformed FC windows for each sub: Pearson's r, gLasso full, and partial covariance windows
%Each type of the FC windows from all subjects are concatenated into one file.
%2. extract the features from the full brain correlation matrix to form the
%featureWin and Z standardize each window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
winSize=69;
step=3; % in TR
numVol=884;
numROI=179;
numWinPerSub=272;
numTotWin=numWinPerSub*numSub;

%define data and resultDir
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
dataDir=[analyDir, '/data/645/all_10min/', covType, '/',session,'/'];
resultDir=[analyDir,'fullBrainFC/results/',covType, '/', session,'/'];
figDir=[analyDir, 'fullBrainFC/figs/GSR/']

% rearrange the regions by lobes, after rearrange:lobe1=F(64ROIs);
% lobe2=T(38ROIs); lobe3=p(35ROIs); lobe4=o(31ROIs);
% lobe5=subcortical(11ROIs)
ROIOrder=xlsread('/home/data/Projects/microstate/DPARSF_preprocessed/doc/Crad_179_lob_hier_label_final.xls','original');
newROIOrder=[find(ROIOrder(:,2)==1); find(ROIOrder(:,2)==2); find(ROIOrder(:,2)==3); find(ROIOrder(:,2)==4); find(ROIOrder(:,2)==5)];
winType={'winFullCor'}; % 'winFullCorLasso'

% for m=1:5
%     length(find(ROIOrder(:,2)==m))
% end

for j=1:length(winType)
    win=char(winType{j})
    %     tmp=load([resultDir,'zfeatureFC_', win, '_', session,'.mat'])
    %     winData=tmp.zfeatureWin;
    %         clustIndx=load([resultDir, 'clustIndx_featureFC_', win, '_', session, '.txt']);
    %     numClust=length(unique(clustIndx));
    
    tmp=load([resultDir,'normFeatureFC_', win, '_', session,'.mat'])
    winData=tmp.normFeature;
    clustIndx=load([resultDir, 'clustIndx_normFeatureFC_', win, '_', session, '.txt']);
    numClust=12;
    
    
    numFeature=size(winData, 2);
    
    
    % plot the cluster mean and rearrange the correlation matrix by lobe
    clustMean=zeros(numClust, numFeature);
    for i=1:numClust
        clust=winData(find(clustIndx==i), :);
        clustMean(i,:)=mean(clust);
    end
    
    % plot the meanClust
    figure(j)
    for i=1:numClust
        
        clustFC=squareform(clustMean(i,:));
        
        % reOrder the columns
        for k=1:numROI
            t=newROIOrder(k)
            colReordered(:,k)=clustFC(:,t);
        end
        
        % reOrder the rows
        for k=1:numROI
            t=newROIOrder(k)
            colAndRowReordered(k,:)=colReordered(t,:);
        end
        
        subplot(3,4,i)
        imagesc(colAndRowReordered)
        colorbar
        caxis([-0.6 0.6])
    end
    %saveas(figure(j), [figDir, session, '_', win, '_stateMap.png'])
end





