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
winType={'winFullCor', 'winFullCorLasso'};

% for m=1:5
%     length(find(ROIOrder(:,2)==m))
% end

for j=1:length(winType)
    win=char(winType{j})
    tmp=load([resultDir,'zfeatureFC_', win, '_', session,'.mat'])
    winData=tmp.zfeatureWin;
    clustIndx=load([resultDir, 'clustIndx_featureFC_', win, '_', session, '.txt']);
    numClust=length(unique(clustIndx));
    numFeature=size(winData, 2);
    
    numWinPerStatePerSub=zeros(numClust, numSub);
    % compute the prct win per state and per sub
    for i=1:numSub
        clustIndxEachSub=clustIndx((1+numWinPerSub*(i-1)):i*numWinPerSub);
        for k=1:numClust
            tmp=length(find(clustIndxEachSub==k));
            numWinPerStatePerSub(k,i)=tmp;
            
        end
    end
    prctWinPerStatePerSub=numWinPerStatePerSub/numWinPerSub*100;
    figure(2)
    imagesc(prctWinPerStatePerSub)
    colorbar
    
    numWinPerState=zeros(numClust);
    for k=1:numClust
        tmp=length(find(clustIndx==k));
        prctWinPerState(k)=tmp/numTotWin*100;
    end
    figure(3)
    bar(prctWinPerState)
    
    saveas(figure(2), [figDir, session, '_', win, '_prctWinPerStatePerSub.png'])
    saveas(figure(3), [figDir, session, '_', win, '_prctWinPerState.png'])
    
    % plot the cluster mean and rearrange the correlation matrix by lobe
    %     clustMean=zeros(numClust, numFeature);
    %     for i=1:numClust
    %         clust=winData(find(clustIndx==i), :);
    %         clustMean(i,:)=mean(clust);
    %     end
    %
    %     % plot the meanClust
    %     figure(j)
    %     for i=1:numClust
    %
    %         clustFC=squareform(clustMean(i,:));
    %
    %         % reOrder the columns
    %         for k=1:numROI
    %             t=newROIOrder(k)
    %             colReordered(:,k)=clustFC(:,t);
    %         end
    %
    %         % reOrder the rows
    %         for k=1:numROI
    %             t=newROIOrder(k)
    %             colAndRowReordered(k,:)=colReordered(t,:);
    %         end
    %
    %         subplot(5,5,i)
    %         imagesc(colAndRowReordered)
    %         colorbar
    %         caxis([-1 1])
    %     end
    %     saveas(figure(j), [figDir, session, '_', win, '_stateMap.png'])
end

% plot the num of cluster per sub
clear
clc
close all
session='session1'
numSub=22
win='FullCor'
numClustPerSub=zeros(numSub,1);
figDir=[analyDir, 'fullBrainFC/figs/GSR/']
clustIndxEachSub=load(['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/results/GSR/', session, '/eachSub/clustIndx_zFeatureWinEachSub_Fullcor.txt']);
for i=1:numSub
    tmp=length(unique(clustIndxEachSub(:, i)));
    numClustPerSub(i,1)=tmp;
end
figure(1)
plot(numClustPerSub)
ylim([12 22])
saveas(figure(1), [figDir, session, '_', win, '_numClustPerSub.png'])



