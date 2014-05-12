
% This script will plot individual surface map and combine individual
% ones into one panel

clear
clc
close all
%%%Get the Surface maps

session='session1';
dataLength='all_10min';
mapType='thresholdedStatsMap';
%mapType='corBrainMapAvg';
surrogateNum=10;

if strcmp(mapType, 'corBrainMapAvg')
    NMin=0; PMin = 0.0000001;
    NMax=-1.5; PMax=2;
else
    NMin=-0.01; PMin = 0.01;
    NMax=-25; PMax=25;
end


Prefix='';

PicturePrefix='';

ClusterSize=0;

SurfaceMapSuffix='_SurfaceMap.jpg';
numClust=5;


ConnectivityCriterion=18;

[BrainNetViewerPath, fileN, extn] = fileparts(which('BrainNet.m'));

SurfFileName=[BrainNetViewerPath,filesep,'Data',filesep,'SurfTemplate',filesep,'BrainMesh_ICBM152_smoothed.nv'];

viewtype='MediumView';

%ColorMap=jet(12);

ColorMap=[1,1,0;1,0.8,0;1,0.6,0;1,0.4118,0;1,0.2667,0;1,0,0;0,0,1;0,0.2667,1;0,0.4118,1;0,0.6,1;0,0.8,1;0,1,1;];
ColorMap=flipdim(ColorMap,1);
cmap1 = colorRamp(ColorMap(1:6,:), 32);
cmap2= colorRamp(ColorMap(7:end,:), 32);
ColorMap=vertcat(cmap1,cmap2)

imgInputDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep,session, '/surrogateData/brainMap/'];

surfaceMapOutputDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/', dataLength, filesep, session,filesep,'surfaceMap', filesep];


for i=1:numClust
    InputName = [imgInputDir,sprintf('thresholdedStatsMap_%dclusters_state%d_surrogate%d_%s_normWin.nii', numClust, i, surrogateNum, session)];
    
    OutputName = [surfaceMapOutputDir,filesep,sprintf('thresholdedStatsMap_5clusters_state%d_surrogate%d_%s_normWin', i, surrogateNum, session),SurfaceMapSuffix];
    
    H_BrainNet = rest_CallBrainNetViewer(InputName,NMin,PMin,ClusterSize,ConnectivityCriterion,SurfFileName,viewtype,ColorMap,NMax,PMax);
    
    eval(['print -r300 -djpeg -noui ''',OutputName,''';']);
end



%%% Auto Draw on a panel

%LeftMaskFile = '/home/data/HeadMotion_YCG/YAN_Scripts/HeadMotion/Parts/Left2FigureMask_BrainNetViewerMediumView.jpg';
close all
DataUpDir = surfaceMapOutputDir;

OutputUpDir = surfaceMapOutputDir;

SurfaceMapSuffix='_SurfaceMap.jpg';


UnitRow = 2168;

%UnitColumn = 1531*2;
UnitColumn = 3095;

BackGroundColor = uint8([255*ones(1,1,3)]);

%imdata_All = 255*ones(UnitRow*length(MeasureSet),UnitColumn*length(ConditionList),3);
cd (DataUpDir)
DirJPG = dir(['*',SurfaceMapSuffix])
totPlots=length(DirJPG);

if totPlots <= 3
    numRow=1;
    numColumn=numImg;
elseif totPlots <=6
    numRow=2;
    numColumn=3;
elseif totPlots <=9
    numRow=3;
    numColumn=3;
else
    numRow=3;
    numColumn=4;
end

imdata_All = repmat(BackGroundColor,[UnitRow*numRow,UnitColumn*numColumn,1]);

%LeftMask = imread(LeftMaskFile);

k=0;
for i=1:numRow
    for  j=1:numColumn
        k=k+1;
        if k<=totPlots
            
            imdata = imread(DirJPG(k).name);
            
            %imdata(LeftMask==255) = 255;
            
            % imdata = imdata(1:end-80,120:1650,:);
            
            imdata=imdata(1:end,120:3214,:);
            
            imdata_All (((i-1)*UnitRow + 1):i*UnitRow,((j-1)*UnitColumn + 1):j*UnitColumn,:) = imdata;
        else
            disp('The total number of plots is smaller than the designed slots')
        end
    end
end

figure

image(imdata_All)

axis off          % Remove axis ticks and numbers

axis image        % Set aspect ratio to obtain square pixels



OutJPGName=[OutputUpDir,filesep,dataLength, '_', session, '_', mapType, '_', seedNum, '_', num2str(totPlots),'clusters.jpg'];

eval(['print -r300 -djpeg -noui ''',OutJPGName,''';']);


