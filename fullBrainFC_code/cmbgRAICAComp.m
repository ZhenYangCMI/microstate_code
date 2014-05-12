
% This script will combine individual surface map for into one panel

clear
clc
close all
%%%Get the Surface maps

clusterMethod='kmeans'
dataDir='/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/testing/figs/modulary/';
OutputUpDir = ['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/testing/figs/'];


numRow=4;
numCol=4;

for i=7:14
    numStates=i;
    
    
    
    % save the img file path
    img=[]
    for j=1:numStates
        img{j,1} = [dataDir, clusterMethod, '_statModulMap_', num2str(numStates), 'clusters_state', num2str(j), '_645TR_0200ROIs_SurfaceMap.jpg'];
    end
    
    % read in the img and save them as one matrixs
    for k=1:numStates
        fileRead=img{k}
        imdata(:, :, :, k)=imread(fileRead);
    end
    
    
    UnitRow = size(imdata, 1); % this num should be corresponding to the size(imdata1, 1)
    
    UnitColumn = size(imdata, 2); % this num should be corresponding to the size(imdata1, 2)
    
    
    BackGroundColor = uint8([255*ones(1,1,3)]);
    imdata_All = repmat(BackGroundColor,[UnitRow*numRow,UnitColumn*numCol,1]);
    
    
    k=0
    for m=1:numRow
        for n=1:numCol
            k=k+1
            if k<numStates+1
                imdata_All (1+(m-1)*UnitRow:m*UnitRow,1+(n-1)*UnitColumn:n*UnitColumn,:) = imdata(:, :, :, k);
            end
        end
    end
    figure
    image(imdata_All)
    axis off          % Remove axis ticks and numbers
    axis image        % Set aspect ratio to obtain square pixels
    OutJPGName=[OutputUpDir,clusterMethod, '_statModulMap_', num2str(numStates), 'states_645TR_0200ROIs.jpg'];
    eval(['print -r300 -djpeg -noui ''',OutJPGName,''';']);
    
end


