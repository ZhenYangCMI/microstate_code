
clear
clc
close all

session='session1';
dataLength='all_10min';
mapType='thresholdedStatsMap';
seedNum='allSeeds';
numStatesSession1=5;
numStatesSession2=6;
numROI=156;

resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,filesep,session,'/'];
maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/smithRSNMask/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep,session, filesep, mapType, filesep, seedNum, '/'];

RSNMaskFile=[maskDir, 'resampled_PNAS_Smith09_rsn10.nii']

[AllVolume, VoxelSize, ImgFileList, Header1, nVolumn] =rest_to4d(RSNMaskFile);
[nDim1 nDim2 nDim3 nRSNs]=size(AllVolume);
brainSize = [nDim1 nDim2 nDim3];

% threshold the RSNs
AllVolume=reshape(AllVolume,[],nRSNs)';
numVoxels=size(AllVolume, 2);
numVoxInRSN=zeros(nRSNs,1);
for i=1:nRSNs
    for j=1:numVoxels
        if AllVolume(i,j)<=3
            AllVolume(i,j)=0;
        else
            AllVolume(i,j)=1;
        end
    end
    numVoxInRSN(i,1)=sum(AllVolume(i,:))
end

% plot the number of voxels in each RSN
% figure(1)
% bar(1:10,numVoxInRSN)
%         title('Number of voxels in a RSN')
%         xlabel('RSN number');
%         ylabel('Number of voxels ');
%      saveas(figure(1), ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep,'NumVoxeInRSN.png'])


% % Convert 2D file back into 4D
% AllVolumeBrain = single(zeros(nRSNs, nDim1*nDim2*nDim3));
% AllVolumeBrain=reshape(AllVolume',[nDim1, nDim2, nDim3, nRSNs]);
%
% % write 4D file as a nift file
% NormAllVolumeBrain=[maskDir,'/','thresholded_Resampled_PNAS_Smith09_rsn10.nii'];
%rest_Write4DNIfTI(AllVolumeBrain,Header1,NormAllVolumeBrain)


if strcmp(session,'session1')
    numStates=numStatesSession1
else
    numStates=numStatesSession2
end
prctOverlapPos=zeros(numStates,nRSNs);
prctOverlapNeg=zeros(numStates,nRSNs);
prctNonOverlap=zeros(numStates,nRSNs);
numVoxOverlapPos=zeros(numStates,nRSNs);
numVoxOverlapNeg=zeros(numStates,nRSNs);
numVoxNonOverlap=zeros(numStates,nRSNs);
numVoxOverlap=zeros(numStates,nRSNs);
k=0;

for m=1:numStates
    for n=1:nRSNs
        
        dataFile=[figDir, 'thresholdedStatsMap_',num2str(numStates),'clusters_state',num2str(m), '_',seedNum,'_',session,'_normWin.nii'];
        [Outdata,VoxDim,Header]=rest_readfile(dataFile);
        [nDim1 nDim2 nDim3]=size(Outdata);
        Outdata=reshape(Outdata,[],1)';
        
        numVoxOverlap(m,n)=length(find(squeeze(AllVolume(n,:))==1&Outdata~=0))
       
        numVoxOverlapPos(m,n)=length(find(squeeze(AllVolume(n,:))==1&Outdata>0))
        numVoxOverlapNeg(m,n)=length(find(squeeze(AllVolume(n,:))==1&Outdata<0))
         if numVoxOverlapPos(m,n)>0
             numVoxOverlapPos(m,n)=numVoxOverlapPos(m,n);
         elseif numVoxOverlapPos(m,n)==0
             numVoxOverlapPos(m,n)=0.00000001;
         else
             disp('Error: numVoxOverlapPos is a negative number.')
         end
             if numVoxOverlapNeg(m,n)>0
             numVoxOverlapNeg(m,n)=numVoxOverlapNeg(m,n);
         elseif numVoxOverlapNeg(m,n)==0
             numVoxOverlapNeg(m,n)=0.00000001;
         else
             disp('Error: numVoxOverlapNeg is a negative number.')
         end
        prctOverlapPos(m,n)=numVoxOverlapPos(m,n)./numVoxInRSN(n)
        prctOverlapNeg(m,n)=numVoxOverlapNeg(m,n)./numVoxInRSN(n)
        
        numVoxNonOverlap(m,n)=numVoxInRSN(n)-numVoxOverlap(m,n)
        prctNonOverlap(m,n)=1-prctOverlapPos(m,n)-prctOverlapNeg(m,n)
        k=k+1
        figure(k)
        x=[prctOverlapPos(m,n),prctOverlapNeg(m,n),prctNonOverlap(m,n)]
        %with the labels for each portion
        %h=pie(x,{'Postivite','Negative','Nonoverlap'});
        % use pieModified the percentage of each slice will not show; use
        % pie function, the percentage label will show
        h=pieModified(x);
        set(h,'LineWidth',3)
        hp = findobj(h, 'Type', 'patch');
        set(hp(1),'FaceColor','r')
        set(hp(2),'FaceColor','b')
        set(hp(3),'FaceColor','w')
        pieFigDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep,session, filesep,'piePlot/individual/'];
        saveas(figure(k), [pieFigDir, session, '_state',num2str(m),'_RSN',num2str(n),'.png'])
        
              
    end
    close all
    
%     plotMatrix=[numVoxNonOverlap(m,:)',numVoxOverlapNeg(m,:)', numVoxOverlapPos(m,:)'];
%     figure(50+m)
%     h=bar(plotMatrix,'stacked');
%     hp = findobj(h, 'Type', 'patch');
%         set(hp(1),'FaceColor','w')
%         set(hp(2),'FaceColor','b')
%         set(hp(3),'FaceColor','r')
%     
%         title(['State',num2str(m),': distribution of overlap and nonoverlap voxels in each RSN'])
%     xlabel('RSN number');
%     ylabel('Number of voxels in a RSN ');
%     set(gca,'YTick',0:1000:8000);
%     ylim([0 8000]);
%     saveas(figure(50+m), ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep,session,filesep,'barPlot',filesep,'state',num2str(m),'_distributOverNonoverlapVoxs.png'])
     
end

% save([resultDir,'prctOverlapPos_',session,'.mat'], 'prctOverlapPos')
% save([resultDir,'numVoxOverlapPos_',session,'.mat'], 'numVoxOverlapPos')
% save([resultDir,'prctOverlapNeg_',session,'.mat'], 'prctOverlapNeg')
% save([resultDir,'numVoxOverlapNeg_',session,'.mat'], 'numVoxOverlapNeg')
close all