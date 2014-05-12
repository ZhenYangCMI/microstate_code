
% This script will plot individual surface map and combine individual
% ones into one panel

clear
clc
close all
%%%Get the Surface maps

session='session2';
dataLength='all_10min';
mapType='thresholdedStatsMap';
%mapType='corBrainMapAvg';
seedList={'seed1','seed2','seed3','seed4'};
numSeed=length(seedList);

for i=1:numSeed
    
seed=char(seedList{i});

if strcmp(mapType, 'corBrainMapAvg')
    NMin=0; PMin = 0.0000001;
    NMax=-1.5; PMax=2;
else
    NMin=-0.00000000001; PMin = 0.00000000001;
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

ColorMap=jet(100);

imgInputDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep,session, filesep, mapType, filesep, seed];

surfaceMapOutputDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/', dataLength, filesep, session,filesep,'surfaceMap', filesep, mapType, filesep, seed];

DirImg = dir([imgInputDir,filesep,Prefix,'*.nii']);
numImg=length(DirImg)

for i=1:numImg
    InputName = [imgInputDir,filesep,DirImg(i).name];
    
    OutputName = [surfaceMapOutputDir,filesep,DirImg(i).name(1:end-4),SurfaceMapSuffix];
    
    H_BrainNet = rest_CallBrainNetViewer(InputName,NMin,PMin,ClusterSize,ConnectivityCriterion,SurfFileName,viewtype,ColorMap,NMax,PMax);
    
    eval(['print -r300 -djpeg -noui ''',OutputName,''';']);
end



%%% Auto Draw on a panel

%LeftMaskFile = '/home/data/HeadMotion_YCG/YAN_Scripts/HeadMotion/Parts/Left2FigureMask_BrainNetViewerMediumView.jpg';
close all
DataUpDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep, session,filesep,'surfaceMap',filesep, mapType, filesep, seed];

OutputUpDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/', dataLength, filesep, session,filesep,'surfaceMap',filesep, 'combined'];

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



OutJPGName=[OutputUpDir,filesep, session, '_', mapType, '_', seed, '_', num2str(totPlots),'clusters.jpg'];

eval(['print -r300 -djpeg -noui ''',OutJPGName,''';']);

end
