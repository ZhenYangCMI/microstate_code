% This script:
% 1.normalizes the time series of each voxel into Z score
% 2.extracts the time series for all seeds and ROIs

clear
clc


subList={'0021002', '0021006', '0021018', '0021024', '1427581', '1793622',...
    '1961098', '2475376', '2799329', '2842950', '3201815', '3313349'...
    '3315657', '3795193', '3808535', '3893245', '4176156', '4288245',...
    '7055197', '8574662', '8735778', '9630905'};
%subList={'0021002'};
sesList={'session1','session2'};
%TRList={'645','2500'};
TRList={'645'};

numSub=size(subList,2);
numSes=size(sesList,2);
numTR=size(TRList,2);
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
maskDir=[analyDir,'mask/'];

for i=1:numTR
    TR=TRList{i};
        TRdir=[analyDir,'data/',TR,'/all_10min'];
    
    for j=1:numSes
        ses=sesList{j};
        if ~exist([TRdir,'/',ses], 'dir')
            mkdir(TRdir, ses)
        end
        sesDir=[TRdir,'/',ses];
        
        for k=1:numSub
            
            disp (['Working on sub ', char(subList{k}),' ......'])
            
            if ~exist([sesDir,'/', char(subList{k})],'dir')
                mkdir(sesDir,char(subList{k}))
            end
            subDir=[sesDir,'/', char(subList{k})];
            
            if strcmp(TR,'645') && strcmp(ses,'session1')
                datadir=['/home/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/preprocessed/FunImgRCW/', char(subList{k}),'/'];
            else
                strcmp(TR,'645') && strcmp(ses,'session2')
                datadir=['/home/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/preprocessed/S2_FunImgRCW/', char(subList{k}),'/'];
                        end
            datafile=[datadir, 'wCovRegressed_4DVolume.nii'];
            
%             % convert the 4D image to 4D matrix
%             [AllVolume, VoxelSize, ImgFileList, Header1, nVolumn] =rest_to4d(datafile);
%             [nDim1 nDim2 nDim3 nDimTimePoints]=size(AllVolume);
%             brainSize = [nDim1 nDim2 nDim3];
            nDim1=61; nDim2=73; nDim3=61;
%             % remove the regions outside of the brain and convert data into 2D
            MaskData=rest_loadmask(nDim1, nDim2, nDim3, [maskDir,'BrainMask_05_61x73x61.img']);
            MaskData =logical(MaskData);%Revise the mask to ensure that it contain only 0 and 1
%             AllVolume=reshape(AllVolume,[],nDimTimePoints)';
%             MaskDataOneDim=reshape(MaskData,1,[]);
%             MaskIndex = find(MaskDataOneDim);
%             AllVolume=AllVolume(:,MaskIndex);
%             
%             
%             % Z_norm the time series for each voxel
%             AllVolume = (AllVolume-repmat(mean(AllVolume),size(AllVolume,1),1))./repmat(std(AllVolume),size(AllVolume,1),1);
%             AllVolume(isnan(AllVolume))=0;
%             
%             
%             % Convert 2D file back into 4D
%             AllVolumeBrain = single(zeros(nDimTimePoints, nDim1*nDim2*nDim3));
%             AllVolumeBrain(:,MaskIndex) = AllVolume;
%             AllVolumeBrain=reshape(AllVolumeBrain',[nDim1, nDim2, nDim3, nDimTimePoints]);
%             
%             
%             % write 4D file as a nift file
             NormAllVolumeBrain=[subDir,'/','norm_AllVolume.nii'];
%             rest_Write4DNIfTI(AllVolumeBrain,Header1,NormAllVolumeBrain)
%             
%             disp ('Time series of each voxel was Z-score normalized.')
            
            
            % extract time series for seeds and ROIs
            
            seed={'ROI Center(mm)=(-4, -45, 8); Radius=3.00 mm.';'ROI Center(mm)=(-3, -45, 23); Radius=3.00 mm.';'ROI Center(mm)=(-2, -39, 27); Radius=3.00 mm.';'ROI Center(mm)=(-2, -36, 35); Radius=3.00 mm.';'ROI Center(mm)=(-3, -51, 65); Radius=3.00 mm.';...
                  'ROI Center(mm)=(-2, -47, 58); Radius=3.00 mm.';'ROI Center(mm)=(-2, -44, 51); Radius=3.00 mm.'; 'ROI Center(mm)=(-2, -31, 45); Radius=3.00 mm.';'ROI Center(mm)=(-2, -60, 59); Radius=3.00 mm.';'ROI Center(mm)=(-2, -56, 51); Radius=3.00 mm.';...
                  'ROI Center(mm)=(-2, -51, 41); Radius=3.00 mm.';'ROI Center(mm)=(-2, -47, 35); Radius=3.00 mm.';'ROI Center(mm)=(-2, -70, 35); Radius=3.00 mm.';'ROI Center(mm)=(-2, -64, 45); Radius=3.00 mm.';'ROI Center(mm)=(-2, -58, 37); Radius=3.00 mm.';...
                  'ROI Center(mm)=(-2, -54, 30); Radius=3.00 mm.';'ROI Center(mm)=(-1, -78, 43); Radius=3.00 mm.'; 'ROI Center(mm)=(-1, -75, 36); Radius=3.00 mm.';'ROI Center(mm)=(-1, -71, 29); Radius=3.00 mm.';'ROI Center(mm)=(-1, -67, 22); Radius=3.00 mm.';...
                  'ROI Center(mm)=(-2, -59, 11); Radius=3.00 mm.'};
            
            [seedROISignals] = y_ExtractROISignal(NormAllVolumeBrain, seed,...
                [subDir,'/21seedROISignal'],MaskData);
            
%             [atlasROISignals] = y_ExtractROISignal(NormAllVolumeBrain, ...
%                 {[maskDir,'final_reduced.nii']},[subDir,'/atlasROISignal'],MaskData,1);
        end
    end
end
disp ('ROI time series extraction done!')

