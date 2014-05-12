
% This script will combine the individual prct overlapping map into one panel

clear
clc
close all
%%%Get the Surface maps

session='session1';
dataLength='all_10min';
mapType='piePlot';

numStatesSession1=5;
numStatesSession2=6;
nRSNs=10;

%%% Auto Draw on a panel

%LeftMaskFile = '/home/data/HeadMotion_YCG/YAN_Scripts/HeadMotion/Parts/Left2FigureMask_BrainNetViewerMediumView.jpg';
close all
DataUpDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep, session,filesep,mapType,filesep, 'individual'];

OutputUpDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/', dataLength, filesep, session,filesep,mapType,filesep, 'combined'];



UnitRow = 821;

UnitColumn = 851;

BackGroundColor = uint8([255*ones(1,1,3)]);


if strcmp(session,'session1')
    numStates=numStatesSession1
else
    numStates=numStatesSession2
end
    
    numRow=numStates;
    numColumn=nRSNs;

imdata_All = repmat(BackGroundColor,[UnitRow*numRow,UnitColumn*numColumn,1]);

%LeftMask = imread(LeftMaskFile);

k=0;
for i=1:numRow
    for  j=1:numColumn
        k=k+1;
                    imgFile=[DataUpDir,filesep,session,'_state',num2str(i),'_RSN',num2str(j),'.png']
            imdata = imread(imgFile);
                                  
            imdata=imdata(40:end-40,200:1050,:);
            
            imdata_All (((i-1)*UnitRow + 1):i*UnitRow,((j-1)*UnitColumn + 1):j*UnitColumn,:) = imdata;
        
        end
    end
totPlots=k;

figure

image(imdata_All)

axis off          % Remove axis ticks and numbers

axis image        % Set aspect ratio to obtain square pixels



OutJPGName=[OutputUpDir,filesep, session, '_', mapType, '_', num2str(k),'subplots.jpg'];

eval(['print -r300 -djpeg -noui ''',OutJPGName,''';']);


