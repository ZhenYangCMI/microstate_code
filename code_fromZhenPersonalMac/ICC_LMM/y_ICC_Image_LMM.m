function y_ICC_Image_LMM(Rate1Dir,Rate2Dir,OutputName,MaskFile)
% function y_ICC_Image_LMM(Rate1Dir,Rate2Dir,OutputName,MaskFile)
% Calculate the Intraclass correlation coefficient for brain images.
% Note: the ICC calculation is based on Xi-Nian Zuo's IPN_icc (http://www.mathworks.com/matlabcentral/fileexchange/22122) which was mainly modified with the Kevin's codes in web. (London kevin.brownhill@kcl.ac.uk)
%   Input:
%     Group1Dir - Cell, directory of the first group. Take average if multiple sessions
%     Group2Dir - Cell, directory of the the second group. Take average if multiple sessions
%   Output:
%     ICC.img - image with Pearson's Correlation Coefficient
%___________________________________________________________________________
% Written by YAN Chao-Gan 110901.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA
% ycg.yan@gmail.com


[ProgramPath, fileN, extn] = fileparts(which('y_ICC_Image_LMM.m'));
addpath([ProgramPath,filesep,'long_mixed_effects_matlab-tools']);
addpath([ProgramPath,filesep,'long_mixed_effects_matlab-tools',filesep,'univariate']);

Rate1Series=0;
for i=1:length(Rate1Dir)
    [Temp,VoxelSize,theImgFileList, Header] =rest_to4d(Rate1Dir{i});
    Rate1Series = Rate1Series + Temp;
    fprintf('\n\tImage Files in Rate %d Directory %d:\n',1, i);
    for itheImgFileList=1:length(theImgFileList)
        fprintf('\t%s%s\n',theImgFileList{itheImgFileList},'.img');
    end
end
Rate1Series = Rate1Series ./ (length(Rate1Dir));


Rate2Series=0;
for i=1:length(Rate2Dir)
    [Temp,VoxelSize,theImgFileList, Header] =rest_to4d(Rate2Dir{i});
    Rate2Series = Rate2Series + Temp;
    fprintf('\n\tImage Files in Rate %d Directory %d:\n',2, i);
    for itheImgFileList=1:length(theImgFileList)
        fprintf('\t%s%s\n',theImgFileList{itheImgFileList},'.img');
    end
end
Rate2Series = Rate2Series ./ (length(Rate2Dir));



[nDim1,nDim2,nDim3,nDim4]=size(Rate2Series);


if ~isempty(MaskFile)
    [MaskData,MaskVox,MaskHead]=rest_readfile(MaskFile);
else
    MaskData=ones(nDim1,nDim2,nDim3);
end

ICCBrain=zeros(nDim1,nDim2,nDim3);


%%%For LMM Only
time = [ones(nDim4,1);2*ones(nDim4,1)];
sID=[[1:nDim4]';[1:nDim4]'];

%%%

for i=1:nDim1
    for j=1:nDim2
        for k=1:nDim3
            if MaskData(i,j,k)
                xA=squeeze(Rate1Series(i,j,k,:));
                xB=squeeze(Rate2Series(i,j,k,:));
                %ICCBrain(i,j,k)=IPN_icc([xA,xB],1,'single');
                
                xAxB = [xA;xB];
                if std(xAxB)~=0  %If no variance, then skip the calculation
                    ICCBrain(i,j,k) = do_ICC([xA;xB], time, [], [], sID);
                else
                    ICCBrain(i,j,k) = 0;
                end
                
            end
        end
    end
end

ICCBrain(isnan(ICCBrain))=0;

rest_WriteNiftiImage(ICCBrain,Header,OutputName);


