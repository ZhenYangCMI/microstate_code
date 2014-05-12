function [KendallWBrain, Header] = y_KendallW_Image(RaterImages, MaskData, AResultFilename)
% Calculate Kendall's W for sets of images. (e.g., different raters or Test Re-Test)
% FORMAT     [KendallWBrain, Header] = y_KendallW_Image(RaterImages, AMaskFilename, AResultFilename)
% Input:
% 	RaterImages     -	Cells of raters (nRater * 1 cells), each rater could be:
%                       1. The directory of 3D image data
%                       2. The filename of one 4D data file
%                       3. a Cell (nFile * 1 cells) of filenames of 3D image data
% 	MaskData		-   the mask file name or 3D mask matrix
%	AResultFilename		the output filename

% Output:
%	KendallWBrain   -   The Kendall's W results
%   Header          -   The NIfTI Header
%	AResultFilename	the filename of Kendall's W result
%-----------------------------------------------------------
% Written by YAN Chao-Gan 120427.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA
% ycg.yan@gmail.com



theElapsedTime =cputime;

fprintf('\n\tKendall''s W computation Start...\n');

[AllVolume,VoxelSize,theImgFileList, Header] =rest_to4d(RaterImages{1});

[nDim1 nDim2 nDim3 nDimTimePoints]=size(AllVolume);
BrainSize = [nDim1 nDim2 nDim3];
VoxelSize = sqrt(sum(Header.mat(1:3,1:3).^2));

if ischar(MaskData)
    if ~isempty(MaskData)
        [MaskData,MaskVox,MaskHead]=rest_readfile(MaskData);
    else
        MaskData=ones(nDim1,nDim2,nDim3);
    end
end

    
MaskDataOneDim=reshape(MaskData,1,[]);
MaskIndex = find(MaskDataOneDim);



%RankSet = zeros(nDimTimePoints,length(MaskIndex),length(RaterImages));
RankSet = repmat(uint16(zeros(nDimTimePoints,1)),[1,length(MaskIndex),length(RaterImages)]);

for iRater = 1:length(RaterImages)

    [AllVolume,VoxelSize,theImgFileList, Header] =rest_to4d(RaterImages{iRater});
    % Convert into 2D
    AllVolume=reshape(AllVolume,[],nDimTimePoints)';
    AllVolume=AllVolume(:,MaskIndex);
    
    
    [AllVolume,SortIndex] = sort(AllVolume,1);
    db=diff(AllVolume,1,1);
    clear AllVolume
    db = db == 0;
    sumdb=sum(db,1);
    
    SortedRanks = repmat(uint16([1:nDimTimePoints]'),[1,length(MaskIndex)]);
    % For those have the same values at the current time point and previous time point (ties)
    if any(sumdb(:))
        TieAdjustIndex=find(sumdb);
        for i=1:length(TieAdjustIndex)
            ranks=SortedRanks(:,TieAdjustIndex(i));
            ties=db(:,TieAdjustIndex(i));
            tieloc = [find(ties); nDimTimePoints+2];
            maxTies = numel(tieloc);
            tiecount = 1;
            while (tiecount < maxTies)
                tiestart = tieloc(tiecount);
                ntied = 2;
                while(tieloc(tiecount+1) == tieloc(tiecount)+1)
                    tiecount = tiecount+1;
                    ntied = ntied+1;
                end
                % Compute mean of tied ranks
                ranks(tiestart:tiestart+ntied-1) = ...
                    sum(ranks(tiestart:tiestart+ntied-1)) / ntied;
                tiecount = tiecount + 1;
            end
            SortedRanks(:,TieAdjustIndex(i))=ranks;
        end
    end
    
    clear db sumdb;
    SortIndexBase = repmat([0:length(MaskIndex)-1].*nDimTimePoints,[nDimTimePoints,1]);
    SortIndex=SortIndex+SortIndexBase;
    clear SortIndexBase;
    I(SortIndex)=SortedRanks;
    clear SortIndex SortedRanks;
    I=uint16(I); 
    
    I = reshape(I,nDimTimePoints,[]);
    
    RankSet(:,:,iRater) = I;
    
end


KendallW = zeros(length(MaskIndex),1);

for iVoxel=1:length(MaskIndex)
    KendallW(iVoxel) = f_kendall(squeeze(RankSet(:,iVoxel,:))); 
end

KendallWBrain=zeros(size(MaskDataOneDim));
KendallWBrain(1,MaskIndex)=KendallW;
KendallWBrain=reshape(KendallWBrain,nDim1, nDim2, nDim3);

Header.pinfo = [1;0;0];
Header.dt    =[16,0];

rest_writefile(KendallWBrain,AResultFilename,BrainSize,VoxelSize,Header,'single');


theElapsedTime =cputime - theElapsedTime;
fprintf('\n\tKendall''s W computation over, elapsed time: %g seconds\n', theElapsedTime);




% calculate kcc for a time series
%---------------------------------------------------------------------------
function B = f_kendall(A)
nk = size(A); n = nk(1); k = nk(2);
SR = sum(A,2); SRBAR = mean(SR);
S = sum(SR.^2) - n*SRBAR^2;
B = 12*S/k^2/(n^3-n);
