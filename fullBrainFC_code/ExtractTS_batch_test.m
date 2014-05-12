% This script:
% 1.normalizes the time series of each voxel into Z score
% 2.extracts the time series for all seeds and ROIs

clear
clc


subList={'3808535','8574662'};
session='session1';


covType='GSR';

numSub=size(subList,2);

analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
maskDir=[analyDir,'fullBrainFC/mask/'];
MaskData=[maskDir,'BrainMask_05_61x73x61.img']

ROIList={'k0050', 'k0100', 'k0200'};

TCDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/testing/TR2500_5min/']


for i=1:length(ROIList)
    numROI=char(ROIList{i})
    for k=1:numSub
        sub=char(subList{k})
        dataDir=[analyDir, '/data/645/all_10min/', covType, '/', session, '/', sub]
        disp (['Working on sub ', char(subList{k}),' ......'])
        
        if ~exist([TCDir,'/', char(subList{k})],'dir')
            mkdir(TCDir,char(subList{k}))
        end
        subDir=[TCDir,'/', char(subList{k})];

datafile=[subDir, '/wCovRegressed_4DVolume.nii'];
            
            % convert the 4D image to 4D matrix
            [AllVolume, VoxelSize, ImgFileList, Header1, nVolumn] =rest_to4d(datafile);
            [nDim1 nDim2 nDim3 nDimTimePoints]=size(AllVolume);
            brainSize = [nDim1 nDim2 nDim3];
            
            % remove the regions outside of the brain and convert data into 2D
            MaskData=rest_loadmask(nDim1, nDim2, nDim3, [maskDir,'BrainMask_05_61x73x61.img']);
            MaskData =logical(MaskData);%Revise the mask to ensure that it contain only 0 and 1
            AllVolume=reshape(AllVolume,[],nDimTimePoints)';
            MaskDataOneDim=reshape(MaskData,1,[]);
            MaskIndex = find(MaskDataOneDim);
            AllVolume=AllVolume(:,MaskIndex);
            
            
            % Z_norm the time series for each voxel
            AllVolume = (AllVolume-repmat(mean(AllVolume),size(AllVolume,1),1))./repmat(std(AllVolume),size(AllVolume,1),1);
            AllVolume(isnan(AllVolume))=0;
            
            
            % Convert 2D file back into 4D
            AllVolumeBrain = single(zeros(nDimTimePoints, nDim1*nDim2*nDim3));
            AllVolumeBrain(:,MaskIndex) = AllVolume;
            AllVolumeBrain=reshape(AllVolumeBrain',[nDim1, nDim2, nDim3, nDimTimePoints]);
            
            
            % write 4D file as a nift file
            NormAllVolumeBrain=[subDir,'/','norm_AllVolume.nii'];
            rest_Write4DNIfTI(AllVolumeBrain,Header1,NormAllVolumeBrain)
            
            disp ('Time series of each voxel was Z-score normalized.')

        
        %NormAllVolumeBrain=[dataDir,'/','norm_AllVolume.nii'];
        
        % extract time series for seeds and ROIs
        output=[numROI, 'signals'];
        
        [output] = y_ExtractROISignal(NormAllVolumeBrain, ...
            {[maskDir, numROI, '.nii']},[subDir,'/TR2500', numROI],MaskData,1);
    end
end

disp ('ROI time series extraction done!')

