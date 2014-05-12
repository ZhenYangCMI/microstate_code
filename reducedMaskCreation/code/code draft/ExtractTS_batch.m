% This script:
% 1.normalizes the time series of each voxel into Z score
% 2.extracts the time series for all seeds and ROIs

clear
clc


SubList={'0021002', '0021006', '0021018', '0021024', '1427581', '1793622',...
    '1961098', '2475376', '2799329', '2842950', '3201815', '3313349'...
    '3315657', '3795193', '3808535', '3893245', '4176156', '4288245',...
    '7055197', '8574662', '8735778', '9630905'};
SesList={'session1','session2'};
TRList={'645','2500'};

N_sub=size(SubList,2);
N_ses=size(SesList,2);
N_TR=size(TRList,2);

maskdir=['/home/data/Projects/microstate/analysis_matlab/mask/'];

for i=1:N_TR
    TR=TRList{i}
    mkdir('/home/data/Projects/microstate/analysis_matlab/data', TR)
    TRdir=['/home/data/Projects/microstate/analysis_matlab/data/',TR];
    for j=1:N_ses
        ses=SesList{j};
        mkdir(TRdir, ses)
        sesdir=[TRdir,'/',ses];
        for k=1:N_sub
            
            disp (['Working on sub ', char(SubList{k}),' ......'])
            
            mkdir(sesdir,char(SubList{k}))
            subdir=[sesdir,'/', char(SubList{k})];
            
            if strcmp(TR,'645') && strcmp(ses,'session1')
                datadir=['/home/data/HeadMotion_YCG/YAN_Work/HeadMotion_YCG/NKI_TRT/RfMRI_mx_645/RegressOutHeadMotion_MNI/HeadMotionRegression_2Friston_24_COV_Global_WM_CSFW/', char(SubList{k}),'/'];
            elseif strcmp(TR,'645'&& strcmp(ses,'session2')
                datadir=['/home/data/HeadMotion_YCG/YAN_Work/HeadMotion_YCG/NKI_TRT/RfMRI_mx_645/S2_RegressOutHeadMotion_MNI/HeadMotionRegression_2Friston_24_COV_Global_WM_CSFW/', char(SubList{k}),'/'];
            elseif strcmp(TR,'2500') && strcmp(ses,'session1')
                datadir=['/home/data/HeadMotion_YCG/YAN_Work/HeadMotion_YCG/NKI_TRT/RfMRI_std_2500/RegressOutHeadMotion_MNI/HeadMotionRegression_2Friston_24_COV_Global_WM_CSFW/', char(SubList{k}),'/'];
            else
                datadir=['/home/data/HeadMotion_YCG/YAN_Work/HeadMotion_YCG/NKI_TRT/RfMRI_std_2500/S2_RegressOutHeadMotion_MNI/HeadMotionRegression_2Friston_24_COV_Global_WM_CSFW/', char(SubList{k}),'/'];
            end
            datafile=[datadir, 'wCovRegressed_4DVolume.nii'];
            
            % convert the 4D image to 4D matrix
            [AllVolume, VoxelSize, ImgFileList, Header1, nVolumn] =rest_to4d(datafile);
            [nDim1 nDim2 nDim3 nDimTimePoints]=size(AllVolume);
            BrainSize = [nDim1 nDim2 nDim3];
            
            % remove the regions outside of the brain and convert data into 2D
            MaskData=rest_loadmask(nDim1, nDim2, nDim3, [maskdir,'BrainMask_05_61x73x61.img']);
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
            NormAllVolumeBrain=[subdir,'/','norm_AllVolume.nii'];
            rest_Write4DNIfTI(AllVolumeBrain,Header1,NormAllVolumeBrain)
            
            disp ('Time series of each voxel was Z-score normalized.')
            
            
            % extract time series for seeds and ROIs
            
            def_seed={'ROI Center(mm)=(-2, -36, 35); Radius=3.00 mm.';'ROI Center(mm)=(-2, -47, 58); Radius=3.00 mm.';...
                'ROI Center(mm)=(-2, -64, 45); Radius=3.00 mm.';'ROI Center(mm)=(-1, -78, 43); Radius=3.00 mm.'};
            
            
            [seed_ROISignals] = y_ExtractROISignal(NormAllVolumeBrain, def_seed,...
                [subdir,'/seed_ROISignal'],MaskData);
            
            [atlas_ROISignals] = y_ExtractROISignal(NormAllVolumeBrain, ...
                {[maskdir,'final_reduced.nii']},[subdir,'/atlas_ROISignal'],MaskData,1);
        end
    end
end
        
        disp ('ROI time series extraction done!')
        
