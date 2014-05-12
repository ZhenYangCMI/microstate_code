function [PCs] = y_Compcor_PC(ADataDir,Nuisance_MaskFilename, OutputName, PCNum, IsNeedDetrend, Band, TR)
% FORMAT [PCs] = y_Compcor_PC(ADataDir,Nuisance_MaskFilename, OutputName, PCNum, IsNeedDetrend, Band, TR)
% Input:
%   ADataDir    -  The data direcotry
%   Nuisance_MaskFilename   -  The Mask file for nuisance area, e.g., the combined mask of WM and CSF
%	OutputName  	-	Output filename
%   PCNum - The number of PCs to be output
%   IsNeedDetrend   -   0: Dot not detrend; 1: Use Matlab's detrend
%   Band            -   Temporal filter band: matlab's ideal filter e.g. [0.01 0.08]
%   TR              -   The TR of scanning. (Used for filtering.)
% Output:
%   PCs - The PCs of the nuisance area (e.g., the combined mask of WM and CSF) for CompCor correction
%__________________________________________________________________________
% Written by YAN Chao-Gan (ycg.yan@gmail.com) on 130808.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA

if ~exist('CUTNUMBER','var')
    CUTNUMBER = 10;
end

if ~exist('PCNum','var')
    PCNum = 5;
end


fprintf('\nExtracting principle components for CompCor Correction:\t"%s"', ADataDir);
[AllVolume,VoxelSize,theImgFileList, Header] =rest_to4d(ADataDir);
[nDim1 nDim2 nDim3 nDimTimePoints]=size(AllVolume);
BrainSize = [nDim1 nDim2 nDim3];

AllVolume=reshape(AllVolume,[],nDimTimePoints)';

if exist('Nuisance_MaskFilename','var') && ~isempty(Nuisance_MaskFilename)
    [MaskData,MaskVox,MaskHead]=rest_readfile(Nuisance_MaskFilename);
else
    MaskData=ones(nDim1,nDim2,nDim3);
end
MaskDataOneDim=reshape(MaskData,1,[]);

AllVolume=AllVolume(:,find(MaskDataOneDim));


% Detrend
if exist('IsNeedDetrend','var') && IsNeedDetrend==1
    %AllVolume=detrend(AllVolume);
    fprintf('\n\t Detrending...');
    SegmentLength = ceil(size(AllVolume,2) / CUTNUMBER);
    for iCut=1:CUTNUMBER
        if iCut~=CUTNUMBER
            Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
        else
            Segment = (iCut-1)*SegmentLength+1 : size(AllVolume,2);
        end
        AllVolume(:,Segment) = detrend(AllVolume(:,Segment));
        fprintf('.');
    end
end

% Filtering
if exist('Band','var') && ~isempty(Band)
    fprintf('\n\t Filtering...');
    SegmentLength = ceil(size(AllVolume,2) / CUTNUMBER);
    for iCut=1:CUTNUMBER
        if iCut~=CUTNUMBER
            Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
        else
            Segment = (iCut-1)*SegmentLength+1 : size(AllVolume,2);
        end
        AllVolume(:,Segment) = y_IdealFilter(AllVolume(:,Segment), TR, Band);
        fprintf('.');
    end
end


% This is for previous meadian angle correction program.
% % Zero temporal Mean and Unit NORM %use std/sqrt(N)  
% AllVolume = (AllVolume-repmat(mean(AllVolume),size(AllVolume,1),1))./repmat(std(AllVolume)*sqrt(nDimTimePoints-1),size(AllVolume,1),1);   %Zero mean and one std
% 
% AllVolume(isnan(AllVolume))=0; %YAN 110123. Set NaN to 0
% This is for previous meadian angle correction program.


% SVD
[U S V] = svd(AllVolume,'econ');

PCs = U(:,1:PCNum);


%Save the results
[pathstr, name, ext] = fileparts(OutputName);

PCs = double(PCs);

save([fullfile(pathstr,[name]), '.mat'], 'PCs')
save([fullfile(pathstr,[name]), '.txt'], 'PCs', '-ASCII', '-DOUBLE','-TABS')

fprintf('\nFinished Extracting principle components for CompCor Correction.\n');





