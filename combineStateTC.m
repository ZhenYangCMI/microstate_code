
% This script will combine the individual prct overlapping map into one panel

clear
clc
close all
%%%Get the Surface maps


dataLength='all_10min';
mapType='stateTC';
sessionList={'session1','session2'};
numSession=length(sessionList);
plotTypeList={'seedsTogether','seedSep'};
numPlotType=length(plotTypeList);
numStatesSession1=5;
numStatesSession2=6;
nRSNs=10;
numSub=22;

%%% Auto Draw on a panel

%LeftMaskFile = '/home/data/HeadMotion_YCG/YAN_Scripts/HeadMotion/Parts/Left2FigureMask_BrainNetViewerMediumView.jpg';
close all
OutputUpDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/', dataLength, filesep, mapType,filesep];



UnitRow = 900;

UnitColumn = 1051;

BackGroundColor = uint8([255*ones(1,1,3)]);




numRow=numSession;
numColumn=numPlotType;

imdata_All = repmat(BackGroundColor,[UnitRow*numRow,UnitColumn*numColumn,1]);

%LeftMask = imread(LeftMaskFile);


for k=1:numSub
    for i=1:numRow
        session=char(sessionList{i});
        DataUpDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep, session,filesep,mapType,filesep, 'individual'];
        
        
        for  j=1:numColumn
            plotType=char(plotTypeList{j});
            
            imgFile=[DataUpDir,filesep,'stateTC_sub',num2str(k),'_', session, plotType,'.png'];
            imdata = imread(imgFile);
            
            imdata=imdata(1:end,50:1100,:);
            
            imdata_All (((i-1)*UnitRow + 1):i*UnitRow,((j-1)*UnitColumn + 1):j*UnitColumn,:) = imdata;
            
        end
    end


    figure
    
    image(imdata_All)
    
    axis off          % Remove axis ticks and numbers
    
    axis image        % Set aspect ratio to obtain square pixels
    
    
    
    OutJPGName=[OutputUpDir,filesep, mapType, '_sub', num2str(k),'combined.jpg'];
    
    eval(['print -r300 -djpeg -noui ''',OutJPGName,''';']);
    
end
close all

