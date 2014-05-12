
% This script will plot individual surface map and combine individual
% ones into one panel

clear
clc
close all
%%%Get the Surface maps

winSize=34;
session='session1';
dataLength='all_10min';
mapType='thresholdedStatsMap';

if (winSize==34)
    numClust=6;
elseif (winSize==136)
    numClust=5;
end


if strcmp(mapType, 'corBrainMapAvg')
    NMin=0; PMin = 0.0000001;
    NMax=-1.5; PMax=2;
else
    NMin=-0.01; PMin = 0.01;
    NMax=-325; PMax=325;
end


Prefix='';

PicturePrefix='';

ClusterSize=0;

SurfaceMapSuffix='_SurfaceMap.jpg';


ConnectivityCriterion=18;

[BrainNetViewerPath, fileN, extn] = fileparts(which('BrainNet.m'));

SurfFileName=[BrainNetViewerPath,filesep,'Data',filesep,'SurfTemplate',filesep,'BrainMesh_ICBM152_smoothed.nv'];

viewtype='MediumView';

%ColorMap=jet(12);

ColorMap=[1,1,0;1,0.9,0;1,0.8,0;1,0.7,0;1,0.6,0;1,0.5,0;0,0.5,1;0,0.6,1;0,0.7,1;0,0.8,1;0,0.9,1;0,1,1;];

ColorMap=flipdim(ColorMap,1);
cmap1 = colorRamp(ColorMap(1:6,:), 32);
cmap2= colorRamp(ColorMap(7:end,:), 32);
ColorMap=vertcat(cmap1,cmap2)

imgInputDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep,session, filesep, mapType, filesep];

surfaceMapOutputDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/', dataLength, filesep, session,'/surfaceMap/difWinSize/'];



for i=1:numClust
    InputName = [imgInputDir,sprintf('thresholdedStatsMap_%dclusters_state%d_allSeeds_%s_normWin_win%d.nii',numClust, i, session, winSize)];
    
    OutputName = [surfaceMapOutputDir,filesep,sprintf('thresholdedStatsMap_%dclusters_state%d_allSeeds_%s_normWin_win%d',numClust, i, session, winSize),SurfaceMapSuffix];
    
    H_BrainNet = rest_CallBrainNetViewer(InputName,NMin,PMin,ClusterSize,ConnectivityCriterion,SurfFileName,viewtype,ColorMap,NMax,PMax);
    
    eval(['print -r300 -djpeg -noui ''',OutputName,''';']);
end



