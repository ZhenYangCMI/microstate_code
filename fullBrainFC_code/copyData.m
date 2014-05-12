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
fileList={'RfMRI_mx_645', 'RfMRI_std_2500'};

for k=1:length(fileList)
    file=char(fileList{k})
    
    for j=1:length(sesList)
        session=char(sesList{j})
        
        for i=1:length(subList)
            sub=char(subList{i});
            
            disp (['Working on sub ', char(subList{i}),' ......'])
            
            funData=sprintf('/home/data/Originals/NKITRT/%s/%s/%s/rest.nii.gz', sub, session, file);
            mkdir(['/home2/data/Projects/microstate/TestGIFT/', file, '/', session, '/FunImg/', sub])
            funDestinationData=['/home2/data/Projects/microstate/TestGIFT/', file, '/', session, '/FunImg/', sub, '/BOLDrest.nii.gz'];
            T1Data=sprintf('/home/data/Originals/NKITRT/%s/anat/mprage.nii.gz', sub)
            mkdir(['/home2/data/Projects/microstate/TestGIFT/', file, '/', session, '/T1Img/', sub])
            T1DestinationDir=['/home2/data/Projects/microstate/TestGIFT/', file, '/', session, '/T1Img/', sub, '/'];
            
            copyfile(funData,funDestinationData)
            copyfile(T1Data,T1DestinationDir)
            
        end
        
    end
    
end