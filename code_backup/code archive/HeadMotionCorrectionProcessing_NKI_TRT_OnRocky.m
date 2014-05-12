% 120912
% Processing on NKI TRT dataset based on preprocessed on IPN server

% Initiate the settings.
% 1. Project Dir
ProjectDir ='/home/data/HeadMotion_YCG/YAN_Work/HeadMotion_YCG/NKI_TRT/RfMRI_mx_645';%
AutoDataProcessParameter.DataProcessDir=[ProjectDir,filesep,'Processing'];
% 2. Subject ID
load('/home/data/HeadMotion_YCG/YAN_Work/HeadMotion_YCG/NKI_TRT/NKITRT_SubID.mat');%
AutoDataProcessParameter.SubjectID=SubID;
AutoDataProcessParameter.SubjectNum=length(AutoDataProcessParameter.SubjectID);
% 3. DARTEL template stored ID
TemplateDir_SubID='0021002';%

% 4. Set Path
addpath /home/data/HeadMotion_YCG/YAN_Program/TRT

[ProgramPath, fileN, extn] = fileparts(which('DPARSFA_run.m'));
Error=[];
addpath([ProgramPath,filesep,'Subfunctions']);
[SPMversion,c]=spm('Ver');
SPMversion=str2double(SPMversion(end));

% 5. Set Common used variable
BrainMaskFile_MNISpace_REST = '/home/data/HeadMotion_YCG/YAN_Program/REST_V1.7_120101/mask/BrainMask_05_61x73x61.img';
load('/home/data/HeadMotion_YCG/YAN_Work/HeadMotion_YCG/Dosenbach_Science_160ROIs_Center.mat');






%%%%%%%
% Then use AllMeasures.m  and  AllMeasuresOnScrubbed.m


rThreshold = 0.15; % p = 0.001 for n = 450 (df = 448, two-tailed) %FOR DEGREE CENTRALITY!!!

TR = 0.645;
IsNeedDetrend = 1;
Filter_Band = [0.01 0.1];
ALFF_Band = [0.01 0.1];




AutoDataProcessParameter.FunctionalSessionNumber=2;
% Multiple Sessions Processing 
% YAN Chao-Gan, 111215 added.
FunSessionPrefixSet={''}; %The first session doesn't need a prefix. From the second session, need a prefix such as 'S2_';
for iFunSession=2:AutoDataProcessParameter.FunctionalSessionNumber
    FunSessionPrefixSet=[FunSessionPrefixSet;{['S',num2str(iFunSession),'_']}];
end





% Get the ID
cd /home/data/Originals/NKITRT
Dir = dir;
for i=3:length(Dir)
    SubID{i-2,1}=Dir(i).name;
end




% Arrange the files
ToDir = '/home/data/HeadMotion_YCG/YAN_Work/HeadMotion_YCG/NKI_TRT/RfMRI_mx_645';
FromDir = '/home/data/Originals/NKITRT';
for i=1:length(SubID)
    mkdir([ToDir,filesep,'FunImg',filesep,SubID{i}]);
    copyfile([FromDir,filesep,SubID{i},filesep,'session1',filesep,'RfMRI_mx_645',filesep,'*'],[ToDir,filesep,'FunImg',filesep,SubID{i}]);
    mkdir([ToDir,filesep,'S2_FunImg',filesep,SubID{i}]);
    copyfile([FromDir,filesep,SubID{i},filesep,'session2',filesep,'RfMRI_mx_645',filesep,'*'],[ToDir,filesep,'S2_FunImg',filesep,SubID{i}]);
    
    mkdir([ToDir,filesep,'T1Img',filesep,SubID{i}]);
    copyfile([FromDir,filesep,SubID{i},filesep,'anat',filesep,'*'],[ToDir,filesep,'T1Img',filesep,SubID{i}]);
    
end








%%%%%%%%%
AutoDataProcessParameter.StartingDirName='FunImg';

AutoDataProcessParameter.FunctionalSessionNumber=2;

% Multiple Sessions Processing 
% YAN Chao-Gan, 111215 added.
FunSessionPrefixSet={''}; %The first session doesn't need a prefix. From the second session, need a prefix such as 'S2_';
for iFunSession=2:AutoDataProcessParameter.FunctionalSessionNumber
    FunSessionPrefixSet=[FunSessionPrefixSet;{['S',num2str(iFunSession),'_']}];
end


AutoDataProcessParameter.RemoveFirstTimePoints=15;
%Remove First Time Points
if (AutoDataProcessParameter.RemoveFirstTimePoints>0)
    for iFunSession=1:AutoDataProcessParameter.FunctionalSessionNumber
        cd([AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName]);
        parfor i=1:AutoDataProcessParameter.SubjectNum
            cd(AutoDataProcessParameter.SubjectID{i});
            
            DirImg=dir('*.nii.gz');  % Search .nii.gz and unzip; YAN Chao-Gan, 120806.
            if length(DirImg)==1
                gunzip(DirImg(1).name);
                delete(DirImg(1).name);
            end
            
            DirImg=dir('*.nii');
            Nii  = nifti(DirImg(1).name);
            
            if size(Nii.dat,4)==900
                y_Write4DNIfTI(Nii.dat(:,:,:,AutoDataProcessParameter.RemoveFirstTimePoints+1:end-1),Nii,DirImg(1).name);
            elseif size(Nii.dat,4)==899
                y_Write4DNIfTI(Nii.dat(:,:,:,AutoDataProcessParameter.RemoveFirstTimePoints+1:end),Nii,DirImg(1).name);
            else
                error(AutoDataProcessParameter.SubjectID{i})
            end
            
            
            cd('..');
            fprintf(['Removing First ',num2str(AutoDataProcessParameter.RemoveFirstTimePoints),' Time Points: ',AutoDataProcessParameter.SubjectID{i},' OK']);
        end
        fprintf('\n');
    end
    %AutoDataProcessParameter.TimePoints=AutoDataProcessParameter.TimePoints-AutoDataProcessParameter.RemoveFirstTimePoints;
end
if ~isempty(Error)
    disp(Error);
    return;
end




%No Slice Timing

AutoDataProcessParameter.StartingDirName='FunImg';
AutoDataProcessParameter.IsRealign=1;
AutoDataProcessParameter.IsCalVoxelSpecificHeadMotion=1;
AutoDataProcessParameter.TimePoints=884;

% Reorient
AutoDataProcessParameter.StartingDirName='FunImgR';
AutoDataProcessParameter.IsNeedReorientT1ImgInteractively=1;
AutoDataProcessParameter.IsNeedReorientFunImgInteractively=1;

%Coregistration and Segment
AutoDataProcessParameter.TimePoints=884;
AutoDataProcessParameter.IsNeedT1CoregisterToFun=1;
AutoDataProcessParameter.IsNeedReorientInteractivelyAfterCoreg=0;
AutoDataProcessParameter.IsSegment=2;
AutoDataProcessParameter.Segment.AffineRegularisationInSegmentation='mni';

%DARTEL
AutoDataProcessParameter.IsDARTEL=1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate FD for each subject

OutputDir=ProjectDir;

for i=1:AutoDataProcessParameter.SubjectNum
    
    DirRP=dir([AutoDataProcessParameter.DataProcessDir,filesep,'RealignParameter',filesep,AutoDataProcessParameter.SubjectID{i},filesep,'rp*']);
    RPFile=[AutoDataProcessParameter.DataProcessDir,filesep,'RealignParameter',filesep,AutoDataProcessParameter.SubjectID{i},filesep,DirRP.name];
    FDThreshold=0.5;
    
    %Calculate FD
    RP=load(RPFile);
    RPDiff=diff(RP);
    RPDiff=[zeros(1,6);RPDiff];
    RPDiffSphere=RPDiff;
    RPDiffSphere(:,4:6)=RPDiffSphere(:,4:6)*50;
    FD=sum(abs(RPDiffSphere),2);
    TemporalMask=ones(length(FD),1);
    Index=find(FD > FDThreshold);
    TemporalMask(Index)=0;
    IndexPrevious=Index;
    for iP=1:1
        IndexPrevious=IndexPrevious-1;
        IndexPrevious=IndexPrevious(IndexPrevious>=1);
        TemporalMask(IndexPrevious)=0;
    end
    IndexNext=Index;
    for iN=1:2
        IndexNext=IndexNext+1;
        IndexNext=IndexNext(IndexNext<=length(FD));
        TemporalMask(IndexNext)=0;
    end
    
    
    FDSet(:,i)=FD;
    FDScrubbedFramesSet_FD05(i)=length(find(TemporalMask==0));
    
end

MeanFD = mean(FDSet)';

FDSet_450=FDSet(1:450,:);
MeanFD_450 = mean(FDSet)';

save([OutputDir,filesep,'FDSet.mat'],'FDSet','MeanFD','FDSet_450','MeanFD_450')





%%%%%%%%%%ALL PROCESSING



%%%%%%%%%
%%%Only use the first 450 time points.
AutoDataProcessParameter.StartingDirName='FunImgR450';

AutoDataProcessParameter.FunctionalSessionNumber=2;

% Multiple Sessions Processing 
% YAN Chao-Gan, 111215 added.
FunSessionPrefixSet={''}; %The first session doesn't need a prefix. From the second session, need a prefix such as 'S2_';
for iFunSession=2:AutoDataProcessParameter.FunctionalSessionNumber
    FunSessionPrefixSet=[FunSessionPrefixSet;{['S',num2str(iFunSession),'_']}];
end


for iFunSession=1:AutoDataProcessParameter.FunctionalSessionNumber
    cd([AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName]);
    parfor i=1:AutoDataProcessParameter.SubjectNum
        cd(AutoDataProcessParameter.SubjectID{i});

        DirImg=dir('rrest.nii');
        Nii  = nifti(DirImg(1).name);
        
        y_Write4DNIfTI(Nii.dat(:,:,:,1:450),Nii,DirImg(1).name);
        
        cd('..');
        
    end
    fprintf('\n');
end





%%%%%%
%Get the normalized and smoothed mean FDvox

parfor i=1:AutoDataProcessParameter.SubjectNum
    for iFunSession=1:AutoDataProcessParameter.FunctionalSessionNumber
        SourceDir=[AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},'VoxelSpecificHeadMotion',filesep,AutoDataProcessParameter.SubjectID{i}];
        
        % Save the mean TDvox and mean FDvox to folder of "MeanVoxelSpecificHeadMotion"
        
        if ~(7==exist([AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'MeanVoxelSpecificHeadMotion_MeanTDvox'],'dir'))
            mkdir([AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'MeanVoxelSpecificHeadMotion_MeanTDvox']);
        end
        movefile([SourceDir,filesep,'MeanTDvox.nii'],[AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'MeanVoxelSpecificHeadMotion_MeanTDvox',filesep,'MeanTDvox_',AutoDataProcessParameter.SubjectID{i},'.nii']);
        
        if ~(7==exist([AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'MeanVoxelSpecificHeadMotion_MeanFDvox'],'dir'))
            mkdir([AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'MeanVoxelSpecificHeadMotion_MeanFDvox']);
        end
        movefile([SourceDir,filesep,'MeanFDvox.nii'],[AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'MeanVoxelSpecificHeadMotion_MeanFDvox',filesep,'MeanFDvox_',AutoDataProcessParameter.SubjectID{i},'.nii']);
    end
end









%Normalize on Results
AutoDataProcessParameter.StartingDirName = 'Results';
AutoDataProcessParameter.IsNormalize=3;
AutoDataProcessParameter.Normalize.Timing='OnResults';
AutoDataProcessParameter.Normalize.BoundingBox=[-90 -126 -72;90 90 108];
AutoDataProcessParameter.Normalize.VoxSize=[3 3 3];

AutoDataProcessParameter.FunctionalSessionNumber=2;

% Multiple Sessions Processing 
% YAN Chao-Gan, 111215 added.
FunSessionPrefixSet={''}; %The first session doesn't need a prefix. From the second session, need a prefix such as 'S2_';
for iFunSession=2:AutoDataProcessParameter.FunctionalSessionNumber
    FunSessionPrefixSet=[FunSessionPrefixSet;{['S',num2str(iFunSession),'_']}];
end


MeasureSet={'MeanVoxelSpecificHeadMotion_MeanTDvox','MeanVoxelSpecificHeadMotion_MeanFDvox'};


fprintf(['Normalizing the resutls into MNI space...\n']);


parfor i=1:AutoDataProcessParameter.SubjectNum
    
    FileList=[];
    for iFunSession=1:AutoDataProcessParameter.FunctionalSessionNumber
        for iMeasure=1:length(MeasureSet)
            cd([AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,filesep,MeasureSet{iMeasure}]);
            DirImg=dir(['*',AutoDataProcessParameter.SubjectID{i},'*.img']);
            for j=1:length(DirImg)
                FileList=[FileList;{[AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,filesep,MeasureSet{iMeasure},filesep,DirImg(j).name]}];
            end
            
            DirImg=dir(['*',AutoDataProcessParameter.SubjectID{i},'*.nii']);
            for j=1:length(DirImg)
                FileList=[FileList;{[AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,filesep,MeasureSet{iMeasure},filesep,DirImg(j).name]}];
            end
        end
        
    end
    
    % Set the mean functional image % YAN Chao-Gan, 120826
    DirMean=dir([AutoDataProcessParameter.DataProcessDir,filesep,'RealignParameter',filesep,AutoDataProcessParameter.SubjectID{i},filesep,'mean*.img']);
    if isempty(DirMean)  %YAN Chao-Gan, 111114. Also support .nii files.
        DirMean=dir([AutoDataProcessParameter.DataProcessDir,filesep,'RealignParameter',filesep,AutoDataProcessParameter.SubjectID{i},filesep,'mean*.nii.gz']);% Search .nii.gz and unzip; YAN Chao-Gan, 120806.
        if length(DirMean)==1
            gunzip([AutoDataProcessParameter.DataProcessDir,filesep,'RealignParameter',filesep,AutoDataProcessParameter.SubjectID{i},filesep,DirMean(1).name]);
            delete([AutoDataProcessParameter.DataProcessDir,filesep,'RealignParameter',filesep,AutoDataProcessParameter.SubjectID{i},filesep,DirMean(1).name]);
        end
        DirMean=dir([AutoDataProcessParameter.DataProcessDir,filesep,'RealignParameter',filesep,AutoDataProcessParameter.SubjectID{i},filesep,'mean*.nii']);
    end
    MeanFilename = [AutoDataProcessParameter.DataProcessDir,filesep,'RealignParameter',filesep,AutoDataProcessParameter.SubjectID{i},filesep,DirMean(1).name];
    
    FileList=[FileList;{MeanFilename}]; %YAN Chao-Gan, 120826. Also normalize the mean functional image.
    
    
    if (AutoDataProcessParameter.IsNormalize==1) %Normalization by using the EPI template directly
        
        SPMJOB = load([ProgramPath,filesep,'Jobmats',filesep,'Normalize.mat']);
        
        SPMJOB.jobs{1,1}.spatial{1,1}.normalise{1,1}.estwrite.subj(1,1).source={MeanFilename};
        SPMJOB.jobs{1,1}.spatial{1,1}.normalise{1,1}.estwrite.subj(1,1).resample=FileList;
        
        [SPMPath, fileN, extn] = fileparts(which('spm.m'));
        SPMJOB.jobs{1,1}.spatial{1,1}.normalise{1,1}.estwrite.eoptions.template={[SPMPath,filesep,'templates',filesep,'EPI.nii,1']};
        SPMJOB.jobs{1,1}.spatial{1,1}.normalise{1,1}.estwrite.roptions.bb=AutoDataProcessParameter.Normalize.BoundingBox;
        SPMJOB.jobs{1,1}.spatial{1,1}.normalise{1,1}.estwrite.roptions.vox=AutoDataProcessParameter.Normalize.VoxSize;
        
        if SPMversion==5
            spm_jobman('run',SPMJOB.jobs);
        elseif SPMversion==8  %YAN Chao-Gan, 090925. SPM8 compatible.
            SPMJOB.jobs = spm_jobman('spm5tospm8',{SPMJOB.jobs});
            spm_jobman('run',SPMJOB.jobs{1});
        else
            uiwait(msgbox('The current SPM version is not supported by DPARSF. Please install SPM5 or SPM8 first.','Invalid SPM Version.'));
            Error=[Error;{['Error in Normalize: The current SPM version is not supported by DPARSF. Please install SPM5 or SPM8 first.']}];
        end
    end
    
    if (AutoDataProcessParameter.IsNormalize==2) %Normalization by using the T1 image segment information
        %Normalize-Write: Using the segment information
        SPMJOB = load([ProgramPath,filesep,'Jobmats',filesep,'Normalize_Write.mat']);
        
        MatFileDir=dir([AutoDataProcessParameter.DataProcessDir,filesep,'T1ImgSegment',filesep,AutoDataProcessParameter.SubjectID{i},filesep,'*seg_sn.mat']);
        MatFilename=[AutoDataProcessParameter.DataProcessDir,filesep,'T1ImgSegment',filesep,AutoDataProcessParameter.SubjectID{i},filesep,MatFileDir(1).name];
        
        SPMJOB.jobs{1,1}.spatial{1,1}.normalise{1,1}.write.subj.matname={MatFilename};
        SPMJOB.jobs{1,1}.spatial{1,1}.normalise{1,1}.write.subj.resample=FileList;
        SPMJOB.jobs{1,1}.spatial{1,1}.normalise{1,1}.write.roptions.bb=AutoDataProcessParameter.Normalize.BoundingBox;
        SPMJOB.jobs{1,1}.spatial{1,1}.normalise{1,1}.write.roptions.vox=AutoDataProcessParameter.Normalize.VoxSize;
        
        if SPMversion==5
            spm_jobman('run',SPMJOB.jobs);
        elseif SPMversion==8  %YAN Chao-Gan, 090925. SPM8 compatible.
            SPMJOB.jobs = spm_jobman('spm5tospm8',{SPMJOB.jobs});
            spm_jobman('run',SPMJOB.jobs{1});
        else
            uiwait(msgbox('The current SPM version is not supported by DPARSF. Please install SPM5 or SPM8 first.','Invalid SPM Version.'));
            Error=[Error;{['Error in Normalize: The current SPM version is not supported by DPARSF. Please install SPM5 or SPM8 first.']}];
        end
        
    end
    
    if (AutoDataProcessParameter.IsNormalize==3) %Normalization by using DARTEL %YAN Chao-Gan, 111111.
        
        SPMJOB = load([ProgramPath,filesep,'Jobmats',filesep,'Dartel_NormaliseToMNI_FewSubjects.mat']);
        
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.fwhm=[0 0 0];
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.preserve=0;
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.bb=AutoDataProcessParameter.Normalize.BoundingBox;
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.vox=AutoDataProcessParameter.Normalize.VoxSize;
        
        DirImg=dir([AutoDataProcessParameter.DataProcessDir,filesep,'T1ImgNewSegment',filesep,AutoDataProcessParameter.SubjectID{1},filesep,'Template_6.*']);
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.template={[AutoDataProcessParameter.DataProcessDir,filesep,'T1ImgNewSegment',filesep,AutoDataProcessParameter.SubjectID{1},filesep,DirImg(1).name]};
        
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.data.subj(1,1).images=FileList;
        
        DirImg=dir([AutoDataProcessParameter.DataProcessDir,filesep,'T1ImgNewSegment',filesep,AutoDataProcessParameter.SubjectID{i},filesep,'u_*']);
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.data.subj(1,1).flowfield={[AutoDataProcessParameter.DataProcessDir,filesep,'T1ImgNewSegment',filesep,AutoDataProcessParameter.SubjectID{i},filesep,DirImg(1).name]};
        
        spm_jobman('run',SPMJOB.matlabbatch);
    end
end


%Copy the Normalized results to ResultsW
for iFunSession=1:AutoDataProcessParameter.FunctionalSessionNumber
    parfor iMeasure=1:length(MeasureSet)
        mkdir([AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,'W',filesep,MeasureSet{iMeasure}])
        movefile([AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,filesep,MeasureSet{iMeasure},filesep,'w*'],[AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,'W',filesep,MeasureSet{iMeasure}])
        fprintf(['Moving Normalized Files:',MeasureSet{iMeasure},' OK']);
    end
    fprintf('\n');
end
AutoDataProcessParameter.StartingDirName=[AutoDataProcessParameter.StartingDirName,'W']; %Now StartingDirName is with new suffix 'W'



%Smooth on Results
AutoDataProcessParameter.IsSmooth=1;
AutoDataProcessParameter.Smooth.Timing='OnResults';
AutoDataProcessParameter.Smooth.FWHM=[4.5 4.5 4.5];

fprintf(['Smoothing the resutls...\n']);



if (AutoDataProcessParameter.IsSmooth==1)
    parfor i=1:AutoDataProcessParameter.SubjectNum
        
        FileList=[];
        for iFunSession=1:AutoDataProcessParameter.FunctionalSessionNumber
            for iMeasure=1:length(MeasureSet)
                cd([AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,filesep,MeasureSet{iMeasure}]);
                DirImg=dir(['*',AutoDataProcessParameter.SubjectID{i},'*.img']);
                for j=1:length(DirImg)
                    FileList=[FileList;{[AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,filesep,MeasureSet{iMeasure},filesep,DirImg(j).name]}];
                end
                
                DirImg=dir(['*',AutoDataProcessParameter.SubjectID{i},'*.nii']);
                for j=1:length(DirImg)
                    FileList=[FileList;{[AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,filesep,MeasureSet{iMeasure},filesep,DirImg(j).name]}];
                end
            end
            
        end
        
        SPMJOB = load([ProgramPath,filesep,'Jobmats',filesep,'Smooth.mat']);
        SPMJOB.jobs{1,1}.spatial{1,1}.smooth.data = FileList;
        SPMJOB.jobs{1,1}.spatial{1,1}.smooth.fwhm = AutoDataProcessParameter.Smooth.FWHM;
        if SPMversion==5
            spm_jobman('run',SPMJOB.jobs);
        elseif SPMversion==8  %YAN Chao-Gan, 090925. SPM8 compatible.
            SPMJOB.jobs = spm_jobman('spm5tospm8',{SPMJOB.jobs});
            spm_jobman('run',SPMJOB.jobs{1});
        else
            uiwait(msgbox('The current SPM version is not supported by DPARSF. Please install SPM5 or SPM8 first.','Invalid SPM Version.'));
        end
        
    end
    
elseif (AutoDataProcessParameter.IsSmooth==2)   %YAN Chao-Gan, 111111. Smooth by DARTEL. The smoothing that is a part of the normalization to MNI space computes these average intensities from the original data, rather than the warped versions. When the data are warped, some voxels will grow and others will shrink. This will change the regional averages, with more weighting towards those voxels that have grows.
    
    parfor i=1:AutoDataProcessParameter.SubjectNum
        
        SPMJOB = load([ProgramPath,filesep,'Jobmats',filesep,'Dartel_NormaliseToMNI_FewSubjects.mat']);
        
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.fwhm=AutoDataProcessParameter.Smooth.FWHM;
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.preserve=0;
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.bb=AutoDataProcessParameter.Normalize.BoundingBox;
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.vox=AutoDataProcessParameter.Normalize.VoxSize;
        
        DirImg=dir([AutoDataProcessParameter.DataProcessDir,filesep,'T1ImgNewSegment',filesep,AutoDataProcessParameter.SubjectID{1},filesep,'Template_6.*']);
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.template={[AutoDataProcessParameter.DataProcessDir,filesep,'T1ImgNewSegment',filesep,AutoDataProcessParameter.SubjectID{1},filesep,DirImg(1).name]};
        
        FileList=[];
        for iFunSession=1:AutoDataProcessParameter.FunctionalSessionNumber
            for iMeasure=1:length(MeasureSet)
                cd([AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName(1:end-1),filesep,MeasureSet{iMeasure}]);
                DirImg=dir(['*',AutoDataProcessParameter.SubjectID{i},'*.img']);
                for j=1:length(DirImg)
                    FileList=[FileList;{[AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName(1:end-1),filesep,MeasureSet{iMeasure},filesep,DirImg(j).name]}];
                end
                
                DirImg=dir(['*',AutoDataProcessParameter.SubjectID{i},'*.nii']);
                for j=1:length(DirImg)
                    FileList=[FileList;{[AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName(1:end-1),filesep,MeasureSet{iMeasure},filesep,DirImg(j).name]}];
                end
            end
            
        end
        
        
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.data.subj(1,1).images=FileList;
        
        DirImg=dir([AutoDataProcessParameter.DataProcessDir,filesep,'T1ImgNewSegment',filesep,AutoDataProcessParameter.SubjectID{i},filesep,'u_*']);
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.data.subj(1,1).flowfield={[AutoDataProcessParameter.DataProcessDir,filesep,'T1ImgNewSegment',filesep,AutoDataProcessParameter.SubjectID{i},filesep,DirImg(1).name]};
        
        spm_jobman('run',SPMJOB.matlabbatch);
        fprintf(['Smooth by using DARTEL:',AutoDataProcessParameter.SubjectID{i},' OK\n']);
    end
    
end


%Copy the Smoothed files to ResultsWS or ResultsS
for iFunSession=1:AutoDataProcessParameter.FunctionalSessionNumber
    parfor iMeasure=1:length(MeasureSet)
        mkdir([AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,'S',filesep,MeasureSet{iMeasure}])
        if (AutoDataProcessParameter.IsSmooth==1)
            movefile([AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,filesep,MeasureSet{iMeasure},filesep,'s*'],[AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,'S',filesep,MeasureSet{iMeasure}])
        elseif (AutoDataProcessParameter.IsSmooth==2) % If smoothed by DARTEL, then the smoothed files still under realign directory.
            movefile([AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName(1:end-1),filesep,MeasureSet{iMeasure},filesep,'s*'],[AutoDataProcessParameter.DataProcessDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,'S',filesep,MeasureSet{iMeasure}])
        end
        fprintf(['Moving Smoothed Files:',MeasureSet{iMeasure},' OK']);
    end
    fprintf('\n');
end

AutoDataProcessParameter.StartingDirName=[AutoDataProcessParameter.StartingDirName,'S']; %Now StartingDirName is with new suffix 'S'














%%%%%%%%%%%%%%%%%%%
%%FOR YANG ZHEN!!!


%%%STEP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regressing out different kind of covariables after head motion regression
% IN ORIGINAL SPACE

DataUpDir=[ProjectDir,filesep,'RegressOutHeadMotion'];
OutputDir=[ProjectDir,filesep,'RegressOutHeadMotion'];

ConditionList={'HeadMotionRegression_2Friston_24'};

for iCondition=1:length(ConditionList)
    
    AutoDataProcessParameter.StartingDirName = ConditionList{iCondition};
    
    
    %2. Regressing out global, WM and CSF
    CovSuffix='_COV_Global_WM_CSF';
    cd([DataUpDir,filesep,AutoDataProcessParameter.StartingDirName]);
    parfor i=1:AutoDataProcessParameter.SubjectNum
        cd(AutoDataProcessParameter.SubjectID{i});
        
        ACovariablesDef=[];
        ACovariablesDef.polort =1;
        ACovariablesDef.ort_file ='';
        ACovariablesDef.CovImgDir={};
        
        BrainMaskFile=[AutoDataProcessParameter.DataProcessDir,filesep,'Masks',filesep,AutoDataProcessParameter.SubjectID{i},'_BrainMask_05_91x109x91.nii'];
        WMMaskFile=[AutoDataProcessParameter.DataProcessDir,filesep,'Masks',filesep,AutoDataProcessParameter.SubjectID{i},'_WhiteMask_09_91x109x91.nii'];
        CSFMaskFile=[AutoDataProcessParameter.DataProcessDir,filesep,'Masks',filesep,AutoDataProcessParameter.SubjectID{i},'_CsfMask_07_91x109x91.nii'];
        GMMaskFile=[AutoDataProcessParameter.DataProcessDir,filesep,'Masks',filesep,AutoDataProcessParameter.SubjectID{i},'_GreyMask_02_91x109x91.nii'];
        
        ACovariablesDef.CovMask = {BrainMaskFile;WMMaskFile;CSFMaskFile};
        
        y_RegressOutImgCovariates([DataUpDir,filesep,AutoDataProcessParameter.StartingDirName,filesep,AutoDataProcessParameter.SubjectID{i}],ACovariablesDef,CovSuffix,'')
        cd ('..');
    end
    
    cd([DataUpDir,filesep,AutoDataProcessParameter.StartingDirName]);
    ToFolderName=[AutoDataProcessParameter.StartingDirName,CovSuffix];
    parfor i=1:AutoDataProcessParameter.SubjectNum
        cd([AutoDataProcessParameter.SubjectID{i}, CovSuffix]);
        mkdir([OutputDir,filesep,ToFolderName,filesep,AutoDataProcessParameter.SubjectID{i}])
        movefile('*',[OutputDir,filesep,ToFolderName,filesep,AutoDataProcessParameter.SubjectID{i}])
        cd('..');
        rmdir([AutoDataProcessParameter.SubjectID{i}, CovSuffix]);
        fprintf(['Moving Coviables Removed Files (',CovSuffix,'):',AutoDataProcessParameter.SubjectID{i},' OK']);
    end
    fprintf('\n');
    
    
    
    
    
    
    
    %3. Regressing out WM and CSF
    CovSuffix='_COV_WM_CSF';
    cd([DataUpDir,filesep,AutoDataProcessParameter.StartingDirName]);
    parfor i=1:AutoDataProcessParameter.SubjectNum
        cd(AutoDataProcessParameter.SubjectID{i});
        
        ACovariablesDef=[];
        ACovariablesDef.polort =1;
        ACovariablesDef.ort_file ='';
        ACovariablesDef.CovImgDir={};
        
        
        BrainMaskFile=[AutoDataProcessParameter.DataProcessDir,filesep,'Masks',filesep,AutoDataProcessParameter.SubjectID{i},'_BrainMask_05_91x109x91.nii'];
        WMMaskFile=[AutoDataProcessParameter.DataProcessDir,filesep,'Masks',filesep,AutoDataProcessParameter.SubjectID{i},'_WhiteMask_09_91x109x91.nii'];
        CSFMaskFile=[AutoDataProcessParameter.DataProcessDir,filesep,'Masks',filesep,AutoDataProcessParameter.SubjectID{i},'_CsfMask_07_91x109x91.nii'];
        GMMaskFile=[AutoDataProcessParameter.DataProcessDir,filesep,'Masks',filesep,AutoDataProcessParameter.SubjectID{i},'_GreyMask_02_91x109x91.nii'];
        
        %ACovariablesDef.CovMask = {BrainMaskFile;WMMaskFile;CSFMaskFile};
        ACovariablesDef.CovMask = {WMMaskFile;CSFMaskFile};
        
        y_RegressOutImgCovariates([DataUpDir,filesep,AutoDataProcessParameter.StartingDirName,filesep,AutoDataProcessParameter.SubjectID{i}],ACovariablesDef,CovSuffix,'')
        cd ('..');
    end
    
    cd([DataUpDir,filesep,AutoDataProcessParameter.StartingDirName]);
    ToFolderName=[AutoDataProcessParameter.StartingDirName,CovSuffix];
    parfor i=1:AutoDataProcessParameter.SubjectNum
        cd([AutoDataProcessParameter.SubjectID{i}, CovSuffix]);
        mkdir([OutputDir,filesep,ToFolderName,filesep,AutoDataProcessParameter.SubjectID{i}])
        movefile('*',[OutputDir,filesep,ToFolderName,filesep,AutoDataProcessParameter.SubjectID{i}])
        cd('..');
        rmdir([AutoDataProcessParameter.SubjectID{i}, CovSuffix]);
        fprintf(['Moving Coviables Removed Files (',CovSuffix,'):',AutoDataProcessParameter.SubjectID{i},' OK']);
    end
    fprintf('\n');
    
    
    
    
end






%STEP
%%%%%
%Normalize 4D Head Motion Corrected Files

ConditionList={'HeadMotionRegression_2Friston_24_COV_Global_WM_CSF','HeadMotionRegression_2Friston_24_COV_WM_CSF'};


ConditionSuffix='';

%Normalize on Results
AutoDataProcessParameter.StartingDirName = 'RegressOutHeadMotion';
AutoDataProcessParameter.IsNormalize=3;
AutoDataProcessParameter.Normalize.Timing='OnFunctionalData';
AutoDataProcessParameter.Normalize.BoundingBox=[-90 -126 -72;90 90 108];
AutoDataProcessParameter.Normalize.VoxSize=[3 3 3];

AutoDataProcessParameter.FunctionalSessionNumber=1;

% Multiple Sessions Processing 
% YAN Chao-Gan, 111215 added.
FunSessionPrefixSet={''}; %The first session doesn't need a prefix. From the second session, need a prefix such as 'S2_';
for iFunSession=2:AutoDataProcessParameter.FunctionalSessionNumber
    FunSessionPrefixSet=[FunSessionPrefixSet;{['S',num2str(iFunSession),'_']}];
end


fprintf(['Normalizing the Data into MNI space...\n']);


parfor i=1:AutoDataProcessParameter.SubjectNum
    
    FileList=[];
    for iFunSession=1:AutoDataProcessParameter.FunctionalSessionNumber
        
        for iCondition=1:length(ConditionList) %%
            cd([ProjectDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,filesep,ConditionList{iCondition},ConditionSuffix,filesep,AutoDataProcessParameter.SubjectID{i}]);
            DirImg=dir(['*.img']);
            for j=1:length(DirImg)
                FileList=[FileList;{[ProjectDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,filesep,ConditionList{iCondition},ConditionSuffix,filesep,AutoDataProcessParameter.SubjectID{i},filesep,DirImg(j).name]}];
            end
            
            DirImg=dir(['*.nii']);
            for j=1:length(DirImg)
                FileList=[FileList;{[ProjectDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,filesep,ConditionList{iCondition},ConditionSuffix,filesep,AutoDataProcessParameter.SubjectID{i},filesep,DirImg(j).name]}];
            end
        end
        
    end
    
    
    if (AutoDataProcessParameter.IsNormalize==3) %Normalization by using DARTEL %YAN Chao-Gan, 111111.
        
        SPMJOB = load([ProgramPath,filesep,'Jobmats',filesep,'Dartel_NormaliseToMNI_FewSubjects.mat']);
        
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.fwhm=[0 0 0];
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.preserve=0;
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.bb=AutoDataProcessParameter.Normalize.BoundingBox;
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.vox=AutoDataProcessParameter.Normalize.VoxSize;
        
        DirImg=dir([AutoDataProcessParameter.DataProcessDir,filesep,'T1ImgNewSegment',filesep,AutoDataProcessParameter.SubjectID{1},filesep,'Template_6.*']);
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.template={[AutoDataProcessParameter.DataProcessDir,filesep,'T1ImgNewSegment',filesep,AutoDataProcessParameter.SubjectID{1},filesep,DirImg(1).name]};
        
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.data.subj(1,1).images=FileList;
        
        DirImg=dir([AutoDataProcessParameter.DataProcessDir,filesep,'T1ImgNewSegment',filesep,AutoDataProcessParameter.SubjectID{i},filesep,'u_*']);
        SPMJOB.matlabbatch{1,1}.spm.tools.dartel.mni_norm.data.subj(1,1).flowfield={[AutoDataProcessParameter.DataProcessDir,filesep,'T1ImgNewSegment',filesep,AutoDataProcessParameter.SubjectID{i},filesep,DirImg(1).name]};
        
        spm_jobman('run',SPMJOB.matlabbatch);
    end
end


%Copy the Normalized results to ResultsW
for iFunSession=1:AutoDataProcessParameter.FunctionalSessionNumber
    for iCondition=1:length(ConditionList) %%
        parfor i=1:AutoDataProcessParameter.SubjectNum
            mkdir([ProjectDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,'_MNI',filesep,ConditionList{iCondition},ConditionSuffix,'W',filesep,AutoDataProcessParameter.SubjectID{i}])
            movefile([ProjectDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,filesep,ConditionList{iCondition},ConditionSuffix,filesep,AutoDataProcessParameter.SubjectID{i},filesep,'w*'],[ProjectDir,filesep,FunSessionPrefixSet{iFunSession},AutoDataProcessParameter.StartingDirName,'_MNI',filesep,ConditionList{iCondition},ConditionSuffix,'W',filesep,AutoDataProcessParameter.SubjectID{i}])
            fprintf(['Moving Normalized Files:',ConditionList{iCondition},' OK']);
        end
    end
    fprintf('\n');
end







%%%STEP
%%%%%%%%%
%Smooth 4D Files after normalization


for iCondition=1:length(ConditionList)
    ConditionList{1,iCondition}=[ConditionList{1,iCondition},'W'];
end

% 
% ConditionList={'FunImgARW','HeadMotionRegression_1Traditional_6W','HeadMotionRegression_2Friston_24W','HeadMotionRegression_3_12W'};
% ConditionList2={'FunImgAR_COV_GlobalW','HeadMotionRegression_1Traditional_6_COV_GlobalW','HeadMotionRegression_2Friston_24_COV_GlobalW','HeadMotionRegression_3_12_COV_GlobalW'};
% ConditionList3={'FunImgAR_COV_Global_WM_CSFW','HeadMotionRegression_1Traditional_6_COV_Global_WM_CSFW','HeadMotionRegression_2Friston_24_COV_Global_WM_CSFW','HeadMotionRegression_3_12_COV_Global_WM_CSFW'};
% 
% ConditionList=[ConditionList,ConditionList2,ConditionList3];
% 


DataUpDir=[ProjectDir,filesep,'RegressOutHeadMotion_MNI'];
OutputDir=[ProjectDir,filesep,'RegressOutHeadMotion_MNI'];
mkdir(OutputDir)



FWHM = [4.5 4.5 4.5];

for iCondition=1:length(ConditionList)
    
    SourceDir_Temp = [DataUpDir,filesep,ConditionList{iCondition}];
    OutpurDir_Temp = [OutputDir,filesep,ConditionList{iCondition},'S'];
    FWHM_Temp = FWHM;
    
    % Too many files, so parfor for each subject
    
    %%%%Smooth Code
    [ProgramPath, fileN, extn] = fileparts(which('DPARSFA_run.m'));
    addpath([ProgramPath,filesep,'Subfunctions']);
    [SPMversion,c]=spm('Ver');
    SPMversion=str2double(SPMversion(end));
    
    parfor i=1:AutoDataProcessParameter.SubjectNum
        SubjectID_Temp = AutoDataProcessParameter.SubjectID{i};
        
        FileList=[];
        DirImg=dir([SourceDir_Temp,filesep,SubjectID_Temp,filesep,'*.img']);
        if isempty(DirImg)
            DirImg=dir([SourceDir_Temp,filesep,SubjectID_Temp,filesep,'*.nii']);
        end
        for j=1:length(DirImg)
            FileList=[FileList;{[SourceDir_Temp,filesep,SubjectID_Temp,filesep,DirImg(j).name]}];
        end
        
        SPMJOB = load([ProgramPath,filesep,'Jobmats',filesep,'Smooth.mat']);
        SPMJOB.jobs{1,1}.spatial{1,1}.smooth.data = FileList;
        SPMJOB.jobs{1,1}.spatial{1,1}.smooth.fwhm = FWHM_Temp;
        if SPMversion==5
            spm_jobman('run',SPMJOB.jobs);
        elseif SPMversion==8  %YAN Chao-Gan, 090925. SPM8 compatible.
            SPMJOB.jobs = spm_jobman('spm5tospm8',{SPMJOB.jobs});
            spm_jobman('run',SPMJOB.jobs{1});
        else
            uiwait(msgbox('The current SPM version is not supported by DPARSF. Please install SPM5 or SPM8 first.','Invalid SPM Version.'));
        end
        
        [SUCCESS,MESSAGE,MESSAGEID] = mkdir([OutpurDir_Temp,filesep,SubjectID_Temp]);
        movefile([SourceDir_Temp,filesep,SubjectID_Temp,filesep,'s*'],[OutpurDir_Temp,filesep,SubjectID_Temp]);
    end
end









