
clear
clc
close all
%%%Get the Surface maps



Prefix='';

PicturePrefix='';



NMin=0;PMin=0.0000001;

NMax=0     ;PMax=6;

ClusterSize=0;

SurfaceMapSuffix='_SurfaceMap.jpg';


ConnectivityCriterion=18;

[BrainNetViewerPath, fileN, extn] = fileparts(which('BrainNet.m'));

SurfFileName=[BrainNetViewerPath,filesep,'Data',filesep,'SurfTemplate',filesep,'BrainMesh_ICBM152_smoothed.nv'];

viewtype='MediumView';

ColorMap=jet(100);

imgInputDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/figs/noGSRFiltered/'];

surfaceMapOutputDir = imgInputDir;

DirImg = dir([imgInputDir,filesep,Prefix,'*.nii']);

for i=1:length(DirImg)
    InputName = [imgInputDir,filesep,DirImg(i).name];
    
    OutputName = [surfaceMapOutputDir,filesep,DirImg(i).name(1:end-4),SurfaceMapSuffix];
    
    H_BrainNet = rest_CallBrainNetViewer(InputName,NMin,PMin,ClusterSize,ConnectivityCriterion,SurfFileName,viewtype,ColorMap,NMax,PMax);
    
    eval(['print -r300 -djpeg -noui ''',OutputName,''';']);
end




