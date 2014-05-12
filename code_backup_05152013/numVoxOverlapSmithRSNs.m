
clear
clc
close all

session='session2';
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
numVoxInRSN=zeros(nRSNs,1)
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


% Convert 2D file back into 4D
AllVolumeBrain = single(zeros(nRSNs, nDim1*nDim2*nDim3));
AllVolumeBrain=reshape(AllVolume',[nDim1, nDim2, nDim3, nRSNs]);

% write 4D file as a nift file
NormAllVolumeBrain=[maskDir,'/','thresholded_Resampled_PNAS_Smith09_rsn10.nii'];
%rest_Write4DNIfTI(AllVolumeBrain,Header1,NormAllVolumeBrain)


if strcmp(session,'session1')
    numStates=numStatesSession1
else
    numStates=numStatesSession2
end
numVoxOverlap=zeros(numStates,nRSNs);
prctOverlap=zeros(numStates,nRSNs);
k=0;

for m=1:numStates
    for n=1:nRSNs
        dataFile=[figDir, 'thresholdedStatsMap_',num2str(numStates),'clusters_state',num2str(m), '_',seedNum,'_',session,'_normWin.nii']
        [Outdata,VoxDim,Header]=rest_readfile(dataFile);
        [nDim1 nDim2 nDim3]=size(Outdata);
        Outdata=reshape(Outdata,[],1)';
        Outdata=logical(Outdata);
        
        numVoxOverlap(m,n)=length(find(squeeze(AllVolume(n,:))==1&Outdata==1))
        prctOverlap(m,n)=length(find(squeeze(AllVolume(n,:))==1&Outdata==1))./numVoxInRSN(n)
        k=k+1;
        %     figure(k)
        %     r=prctOverlap(m,n)*10;
        %     filledCircle([10,10],r,100,'g');
        %     axis([0 20 0 20])
        %     axis off
        %     circleFigDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep,session, filesep,'circleMap/individual/'];
        %saveas(figure(k), [circleFigDir, session, '_state',num2str(m),'_RSN',num2str(n),'.png'])
    end
    close all
    figure(70+m)
    bar(1:10,numVoxOverlap(m,:))
    title(['State',num2str(m),':number of voxels overlap with RSNs'])
    xlabel('RSN number');
    ylabel('Number of overlapping voxels ');
    set(gca,'YTick',0:1000:8000);
    ylim([0 8000]);
    saveas(figure(70+m), ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep,session,'state',num2str(m),'_NumVoxeOverlapping.png'])
    
end

save([resultDir,'prctOverlap_',session,'.mat'], 'prctOverlap')
save([resultDir,'numVoxOverlap_',session,'.mat'], 'numVoxOverlap')
close all