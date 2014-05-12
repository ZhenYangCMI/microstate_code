% This script:
% 1.normalizes the time series of each voxel into Z score
% 2.extracts the time series for all seeds and ROIs

clear
clc


subList={'0021002', '0021006', '0021018', '0021024', '1427581', '1793622',...
    '1961098', '2475376', '2799329', '2842950', '3201815', '3313349'...
    '3315657', '3795193', '3808535', '3893245', '4176156', '4288245',...
    '7055197', '8574662', '8735778', '9630905'};
sesList={'session1','session2'};
%sesList={'session2'};
%TRList={'645','2500'};
TRList={'645'};
covType='GSR';

numSub=size(subList,2);
numSes=size(sesList,2);
numTR=size(TRList,2);
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
maskDir=[analyDir,'fullBrainFC/mask/'];
MaskData=[maskDir,'BrainMask_05_61x73x61.img']

for i=1:numTR
    TR=TRList{i};
        TRdir=[analyDir,'data/',TR,'/all_10min/', covType];
    
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
            
            NormAllVolumeBrain=[subDir,'/','norm_AllVolume.nii'];
            
            % extract time series for seeds and ROIs
            
            [cc179ROISignals] = y_ExtractROISignal(NormAllVolumeBrain, ...
                {[maskDir,'Crad179.nii']},[subDir,'/cc179ROISignals'],MaskData,1);
        end
    end
end
disp ('ROI time series extraction done!')

