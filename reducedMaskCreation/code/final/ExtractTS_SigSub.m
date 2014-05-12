% This script z score normalize the time series for each voxel and extract
% the time series for all seeds and ROIs

clear
clc



% Z normalize the time series of each voxel
[AllVolume, VoxelSize, ImgFileList, Header1, nVolumn] =rest_to4d('/Users/zhenyang/Documents/microstate/data/0021018/swCovRegressed_4DVolume');
[nDim1 nDim2 nDim3 nDimTimePoints]=size(AllVolume);
BrainSize = [nDim1 nDim2 nDim3];

% Convert into 2D

MaskData=rest_loadmask(nDim1, nDim2, nDim3, '/Users/zhenyang/Documents/microstate/mask/mask_used/BrainMask_05_61x73x61.img');

MaskData =logical(MaskData);%Revise the mask to ensure that it contain only 0 and 1


AllVolume=reshape(AllVolume,[],nDimTimePoints)';

MaskDataOneDim=reshape(MaskData,1,[]);
MaskIndex = find(MaskDataOneDim);
 
AllVolume=AllVolume(:,MaskIndex);


% Z_norm the time series for each voxel
AllVolume = (AllVolume-repmat(mean(AllVolume),size(AllVolume,1),1))./repmat(std(AllVolume),size(AllVolume,1),1);   
AllVolume(isnan(AllVolume))=0;


% Convert into 4D
AllVolumeBrain = single(zeros(nDimTimePoints, nDim1*nDim2*nDim3));
AllVolumeBrain(:,MaskIndex) = AllVolume;
 
AllVolumeBrain=reshape(AllVolumeBrain',[nDim1, nDim2, nDim3, nDimTimePoints]);


% write 4D file as a nift file
rest_Write4DNIfTI(AllVolumeBrain,Header1,'/Users/zhenyang/Documents/microstate/data/0021018/norm_AllVolume.nii')

disp ('Time series of each voxel of Z score normalized.')

 
% extract Seed and ROI time series

def_seed={'ROI Center(mm)=(-2, -36, 35); Radius=3.00 mm.';'ROI Center(mm)=(-2, -47, 58); Radius=3.00 mm.';...
'ROI Center(mm)=(-2, -64, 45); Radius=3.00 mm.';'ROI Center(mm)=(-1, -78, 43); Radius=3.00 mm.'}

[seed_ROISignals] = y_ExtractROISignal('/Users/zhenyang/Documents/microstate/data/0021018/norm_AllVolume', def_seed,...
'/Users/zhenyang/Documents/microstate/data/time_series/seed_ROISignal',MaskData);

[atlas_ROISignals] = y_ExtractROISignal('/Users/zhenyang/Documents/microstate/data/0021018/norm_AllVolume', ...
{'/Users/zhenyang/Documents/microstate/mask/final_ROIs_Crad179/final_reduced.nii'},...
'/Users/zhenyang/Documents/microstate/data/time_series/atlas_ROISignal',MaskData,1);

disp ('Extract the time series for the seeds and the selected ROIs.')