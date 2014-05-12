% PCA analysis
clear
clc

% Global variable
session='session1'
maskDir='/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/mask/';
subList={'0021002', '0021006', '0021018', '0021024', '1427581', '1793622',...
    '1961098', '2475376', '2799329', '2842950', '3201815', '3313349'...
    '3315657', '3795193', '3808535', '3893245', '4176156', '4288245',...
    '7055197', '8574662', '8735778', '9630905'};

ROIoutputDir='/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/data/noGSRFiltered/';
resultDir='/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/results/noGSRFiltered/'
figDir = '/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/figs/noGSRFiltered/';

numVol=884;
winSize=69;
step=3;
numWin=floor((numVol-winSize)/step)+1;
numSub=length(subList);
numWinAllSub=numWin*numSub;



%% 1. mask creation AAL mask with ROI 1 to 88

% mask=[mask, '/AAL_61x73x61_YCG.nii']
% [MaskData,MaskVox,MaskHead]=rest_readfile(mask);
%
% MaskData(find(MaskData>88))=0;
% reducedMask=MaskData;
%
% fileName='/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/mask/AAL_61x73x61_YCG_reduced.nii';
% rest_WriteNiftiImage(reducedMask,MaskHead,fileName);


%% 2. Extract the TS from ROIs

for i=1:length(subList)
    sub=char(subList{i})
    
    if ~exist([ROIoutputDir, sub], 'dir')
        mkdir(ROIoutputDir,sub)
    end
    subDir=[ROIoutputDir,sub];
    
    data=['/home/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/preprocessed/noGSR/FunImgRCFW/', sub, '/wFiltered_4DVolume.nii'];
    
    [AllVolume, VoxelSize, ImgFileList, Header1, nVolumn] =rest_to4d(data);
    [nDim1 nDim2 nDim3 nDimTimePoints]=size(AllVolume);
    brainSize = [nDim1 nDim2 nDim3];
    
    % remove the regions outside of the brain and convert data into 2D
    MaskData=rest_loadmask(nDim1, nDim2, nDim3, [maskDir,'BrainMask_05_61x73x61.img']);
    MaskData =logical(MaskData);%Revise the mask to ensure that it contain only 0 and 1
    AllVolume=reshape(AllVolume,[],nDimTimePoints)';
    MaskDataOneDim=reshape(MaskData,1,[]);
    MaskIndex = find(MaskDataOneDim);
    AllVolume=AllVolume(:,MaskIndex);
    
    
    % Z_norm the time series for each voxel
    AllVolume = (AllVolume-repmat(mean(AllVolume),size(AllVolume,1),1))./repmat(std(AllVolume),size(AllVolume,1),1);
    AllVolume(isnan(AllVolume))=0;
    
    
    % Convert 2D file back into 4D
    AllVolumeBrain = single(zeros(nDimTimePoints, nDim1*nDim2*nDim3));
    AllVolumeBrain(:,MaskIndex) = AllVolume;
    AllVolumeBrain=reshape(AllVolumeBrain',[nDim1, nDim2, nDim3, nDimTimePoints]);
    
    
    % write 4D file as a nift file
    NormAllVolumeBrain=[subDir,'/','norm_AllVolume.nii'];
    rest_Write4DNIfTI(AllVolumeBrain,Header1,NormAllVolumeBrain)
    
    disp ('Time series of each voxel was Z-score normalized.')
    
    
    % extract time series for seeds and ROIs
    
    [aalROISignals] = y_ExtractROISignal(NormAllVolumeBrain, ...
        {[maskDir,'AAL_61x73x61_YCG_reduced.nii']},[subDir,'/aalROISignal'],MaskData,1);
end

disp ('ROI time series extraction done!')

%% 3. FC window creation


ROISignalName='ROISignals_aalROISignal.mat'

% create the convolved window
[finalWin]=winCreation(winSize,0);
finalWin=finalWin';

% segment the TS into TS windows

for k=1:numSub
    sub=char(subList{k});
    disp (['Working on sub ', sub,' ......'])
    subDir=[ROIoutputDir,sub];
    ROISignals = load([subDir,'/', ROISignalName]);
    TC=ROISignals.ROISignals;
    
    % apply the sliding window to the time series
    asize = size(TC);
    
    % win(winSize,asize(2),numWinAllSub);
    for q=1+numWin*(k-1):numWin*k;
        for n=((q-1)-(k-1)*numWin)*step+1:winSize+((q-1)-(k-1)*numWin)*step
            for m=1:asize(2)
                win(n-((q-1)-(k-1)*numWin)*step,m,q)=TC(n,m)*finalWin((n-((q-1)-(k-1)*numWin)*step),1);
            end
        end
    end
end
disp (['Window applying done for all subjects. All windows of time series are saved in one matrix!'])

% generate full correlation for each window
disp (['Compute the full correlation for each window'])
for q=1:numWinAllSub
    
    % compute the full correlation: Pearson's r
    fullCorrWin(:,:,q)=corrcoef(win(:,:,q));
end

% Fisher z tranform the correlations
zFullCorrWin=0.5*log((1+fullCorrWin)./(1-fullCorrWin));

disp('Full correlation are computed for each window.')

save([resultDir, 'fullCorrWin_', session, '.mat'], 'fullCorrWin')
save([resultDir, 'zFullCorrWin_', session, '.mat'], 'zFullCorrWin')

%% 4. Extract feature
tmp=load([resultDir, 'zFullCorrWin_', session, '.mat'])
FCWinFullBrain=tmp.zFullCorrWin;
[featureWin, zFeatureWin] = featureExtract(FCWinFullBrain);
save([resultDir, 'feature_', session, '.mat'], 'featureWin')
save([resultDir, 'zNormWinFeature_', session, '.mat'], 'zFeatureWin')

%% 5 standardize the feature
tmp=load([resultDir, '/feature_', session, '.mat'])
feature=tmp.featureWin;
normFeature1=zeros(size(feature,1), size(feature,2));
normFeature2=zeros(size(feature,1), size(feature,2));
for i=1:numSub
    sub=num2str(i)
    featureSub=feature(1+numWin*(i-1):numWin*i, :);
    featureSub1D=reshape(featureSub, [],1);
    featureSubNorm1=(featureSub-repmat(mean(featureSub1D), size(featureSub,1), size(featureSub,2)))./repmat(std(featureSub1D), size(featureSub,1), size(featureSub,2));
    featureSubNorm2=featureSubNorm1-repmat(mean(featureSubNorm1),size(featureSub, 1), 1);
    normFeature1(1+numWin*(i-1):numWin*i, :)=featureSubNorm1;
    normFeature2(1+numWin*(i-1):numWin*i, :)=featureSubNorm2;
end

save([resultDir,'/normFeature1_',session,'.mat'],'normFeature1')
save([resultDir,'/normFeature2_',session,'.mat'],'normFeature2')

%% 6 run PCA

norm='norm1';
if strcmp(norm, 'norm1')
    loadFile=load([resultDir,'/normFeature1_',session,'.mat'])
    normFeature=loadFile.normFeature1;
else
    loadFile=load([resultDir,'/normFeature2_',session,'.mat'])
    normFeature=loadFile.normFeature2;
end

[COEFF,SCORE,latent, tsquare] = princomp(normFeature);
cumVar=cumsum(latent)./sum(latent);
save([resultDir,'/eigenvector_',norm, '_',session,'.mat'],'COEFF')
save([resultDir,'/eigenvalue_',norm, '_',session,'.mat'],'latent')
save([resultDir,'/cumVar_', norm,  '_',session,'.mat'],'cumVar')
save([resultDir,'/tsquare_', norm,  '_',session,'.mat'],'tsquare')

x=1:length(latent);
close all
figure(1)
subplot(2,1,1)
[AX, H1, H2]=plotyy(x,latent, x, cumVar);
set(get(AX(1),'Ylabel'),'String','Eigen Value')
set(get(AX(2),'Ylabel'),'String','Cummulative Variance')
xlabel('Number of Eigen Values')
title('All Number of Eigen Values')
set(H1,'LineStyle','-', 'LineWidth', 2)
set(H2,'LineStyle','-', 'LineWidth', 2)

num=100
subplot(2,1,2)
[AX, H1, H2]=plotyy(x(1:num),latent(1:num), x(1:num), cumVar(1:num));
set(get(AX(1),'Ylabel'),'String','Eigen Value')
set(get(AX(2),'Ylabel'),'String','Cummulative Variance')
xlabel('Number of Eigen Values')
title(['First ', num2str(num),  'Eigen Values'])
set(H1,'LineStyle','-', 'LineWidth', 2)
set(H2,'LineStyle','-', 'LineWidth', 2)
saveas(figure(1), [figDir, 'eigenValue_', norm, '_', session, '.jpg'])

numPC=10
figure(2)
for i=1:numPC
    PC=squareform(COEFF(:,i));
    subplot(3,4,i)
    imagesc(PC)
    colorbar
if i==1
    caxis([-0.025 0.025])
else
caxis([-0.04 0.04])
end
    title(['PC', num2str(i)])
    axis square
end
saveas(figure(2), [figDir, 'PCs_', norm, '_', session '.jpg'])


