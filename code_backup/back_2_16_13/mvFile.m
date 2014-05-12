

dataDir=['/home/data/HeadMotion_YCG/YAN_Work/HeadMotion_YCG/NKI_TRT/RfMRI_mx_645/Processing/T1ImgNewSegment/'];
saveDir=['/home/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/T1Img/'];
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;
numSub=length(subList);

for i=1:numSub
    subDataDir=[dataDir,subList{i},'/'];
    disp (subDataDir)
    mkdir(saveDir,subList{i})
    subSaveDir=[saveDir,subList{i},'/'];
    disp (subSaveDir)
    copyfile([subDataDir,'mprage.nii'],subSaveDir)
end


        