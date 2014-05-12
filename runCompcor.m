clear
clc

subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;
numSub=length(subList);

PCNum=5;
IsNeedDetrend=0;
Band=[0.01 0.1];
TR=0.645;

for i=1:numSub
sub=subList(i)
ADataDir=['/home2/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/preprocessed/S2_FunImgR/', char(sub)];
Nuisance_MaskFilename=['/home2/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/preprocessed/Masks/', char(sub), '_CsfWhiteMask_91x109x91'];
OutputName=['/home2/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/preprocessed/S2_compCorPC/', char(sub), '_', num2str(PCNum), 'componets.txt'];
[PCs] = y_Compcor_PC(ADataDir,Nuisance_MaskFilename, OutputName, PCNum, IsNeedDetrend, Band, TR);
end