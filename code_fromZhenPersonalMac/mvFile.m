
% copy the anatomic data over from Chao-gan's folder

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


%  copy the raw functional data over from Chao-gan's folder
clear
clc

sourceDir=['/home/data/HeadMotion_YCG/YAN_Work/HeadMotion_YCG/NKI_TRT/RfMRI_mx_645/Processing/FunImgR_copy/'];
destinatDir=['/home/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/FunImg'];
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;
numSub=length(subList);


for i=1:numSub
    subSourceDir=[sourceDir,subList{i},'/'];
    disp (subSourceDir)
        mkdir(destinatDir,subList{i})
    subDestinatDir=[destinatDir,'/', subList{i},'/'];
    disp (subDestinatDir)
    copyfile([subSourceDir,'rest.nii'],subDestinatDir)
end


        