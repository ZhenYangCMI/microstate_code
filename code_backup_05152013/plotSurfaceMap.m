
clear
clc
close all
%%%Get the Surface maps



Prefix='';

PicturePrefix='';



NMin=-0.000001;PMin=0.000001;

NMax=-2     ;PMax=2;

ClusterSize=0;

SurfaceMapSuffix='_SurfaceMap.jpg';


ConnectivityCriterion=18;

[BrainNetViewerPath, fileN, extn] = fileparts(which('BrainNet.m'));

SurfFileName=[BrainNetViewerPath,filesep,'Data',filesep,'SurfTemplate',filesep,'BrainMesh_ICBM152_smoothed.nv'];

viewtype='MediumView';

ColorMap=jet(100);


session='session2'
imgInputDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/first_5min/',session,filesep,session,'_ward_euclidean/lambdaOptimalPerSub/corBrainMapAvg/normWin'];

surfaceMapOutputDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/first_5min/',session,filesep,session,'_ward_euclidean/lambdaOptimalPerSub/surfaceMap/normWin'];

DirImg = dir([imgInputDir,filesep,Prefix,'*.nii']);

for i=1:length(DirImg)
    InputName = [imgInputDir,filesep,DirImg(i).name];
    
    OutputName = [surfaceMapOutputDir,filesep,DirImg(i).name(1:end-4),SurfaceMapSuffix];
    
    H_BrainNet = rest_CallBrainNetViewer(InputName,NMin,PMin,ClusterSize,ConnectivityCriterion,SurfFileName,viewtype,ColorMap,NMax,PMax);
    
    eval(['print -r300 -djpeg -noui ''',OutputName,''';']);
end



%%% Auto Draw on a panel

session='session1';
%LeftMaskFile = '/home/data/HeadMotion_YCG/YAN_Scripts/HeadMotion/Parts/Left2FigureMask_BrainNetViewerMediumView.jpg';

DataUpDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/first_5min/',session,filesep,session,'_ward_euclidean/lambdaOptimalPerSub/surfaceMap/normWin'];

OutputUpDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/first_5min/',session,filesep,session,'_ward_euclidean/lambdaOptimalPerSub/surfaceMap/normWin'];

SurfaceMapSuffix='_SurfaceMap.jpg';


UnitRow = 2168;

%UnitColumn = 1531*2;
UnitColumn = 3095;

BackGroundColor = uint8([255*ones(1,1,3)]);

%imdata_All = 255*ones(UnitRow*length(MeasureSet),UnitColumn*length(ConditionList),3);

numRow=2;
numColumn=2;
imdata_All = repmat(BackGroundColor,[UnitRow*numRow,UnitColumn*numColumn,1]);

%LeftMask = imread(LeftMaskFile);

cd (DataUpDir)
DirJPG = dir(['*',SurfaceMapSuffix])
k=0;
totPlots=length(DirJPG);
for i=1:numRow
    for  j=1:numColumn
        k=k+1;
imdata = imread(DirJPG(k).name); 

%imdata(LeftMask==255) = 255;

% imdata = imdata(1:end-80,120:1650,:);

imdata=imdata(1:end,120:3214,:);

imdata_All (((i-1)*UnitRow + 1):i*UnitRow,((j-1)*UnitColumn + 1):j*UnitColumn,:) = imdata;
    end
end

figure

image(imdata_All)

axis off          % Remove axis ticks and numbers

axis image        % Set aspect ratio to obtain square pixels

 

OutJPGName=[OutputUpDir,filesep,num2str(totPlots),'clusters_',session,'_normWin.jpg'];

eval(['print -r300 -djpeg -noui ''',OutJPGName,''';']);


