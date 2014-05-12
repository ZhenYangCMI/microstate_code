
clear
clc
close all


maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/smithRSNMask/'];
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
            AllVolume(i,j)=AllVolume(i,j);
        end
    end
    numVoxInRSN(i,1)=length(find(AllVolume(i,:)~=0))
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
NormAllVolumeBrain=[maskDir,'/','thresholded_Resampled_PNAS_Smith09_rsn10_zscore.nii'];
rest_Write4DNIfTI(AllVolumeBrain,Header1,NormAllVolumeBrain)

