% This script standardize the dyanmic FC matrix for each subject by:
% 1: removing the global mean and dividing by its std
% 2: subtract the row wise mean
clear
clc

session='session1';
winType='winFullCor' % can be winFullCorLasso winFullCor

analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
resultDir=[analyDir,'fullBrainFC/results/GSR/', session];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/figs/GSR/PCA/'];

% plot fig for one sub and compute the mean dyanmic FC

tmp=load([resultDir, '/normFeature2_', winType, '_', session, '.mat'])
feature=tmp.normFeature2;
sub=20; % sub20 with lowest motion
sub20=feature((1+19*272):272*20, :)';
figure(1)
imagesc(sub20)
colorbar
caxis([-1 1])
saveas(figure(1), [figDir, 'sub20_dynamicFC', session '.jpg'])

ROIOrder=xlsread('/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/doc/Cradd179_ROIorder.xlsx','Sheet1','a2:b180');
orderByLobe=ROIOrder(:,1);


i=0
figure(2)
for win=1:30:241
    i=i+1
dataPlot=squareform(sub20(:,win));
tmpReorderCol=dataPlot(:, orderByLobe);
    tmpReorderColAndRow=tmpReorderCol(orderByLobe, :);
    subplot(3,3,i)
imagesc(tmpReorderColAndRow)
colorbar
caxis([-1.5 1.5])
title(['win', num2str(win)])
end
saveas(figure(2), [figDir, 'sub20_norm2dynamicFC_exampleWin_',  session '.jpg'])

% plot average dynamic FC
tmp=load([resultDir, '/featureFC_', winType, '_', session, '.mat'])
feature=tmp.featureWin;
meanFeature=mean(feature);
dataPlot=squareform(meanFeature);
tmpReorderCol=dataPlot(:, orderByLobe);
    tmpReorderColAndRow=tmpReorderCol(orderByLobe, :);
    figure(3)
    imagesc(tmpReorderColAndRow)
colorbar
caxis([-0.3 0.3])
saveas(figure(3), [figDir, 'averageDynamicFC_', session '.jpg'])

%% standardize the feature
tmp=load([resultDir, '/featureFC_', winType, '_', session, '.mat'])
feature=tmp.featureWin;
numWinPerSub=272;
numSub=22
normFeature1=zeros(size(feature,1), size(feature,2));
normFeature2=zeros(size(feature,1), size(feature,2));
for i=1:numSub
    sub=num2str(i)
    featureSub=feature(1+numWinPerSub*(i-1):numWinPerSub*i, :);
    featureSub1D=reshape(featureSub, [],1);
    featureSubNorm1=(featureSub-repmat(mean(featureSub1D), size(featureSub,1), size(featureSub,2)))./repmat(std(featureSub1D), size(featureSub,1), size(featureSub,2));
    featureSubNorm2=featureSubNorm1-repmat(mean(featureSubNorm1),size(featureSub, 1), 1);
    normFeature1(1+numWinPerSub*(i-1):numWinPerSub*i, :)=featureSubNorm1;
    normFeature2(1+numWinPerSub*(i-1):numWinPerSub*i, :)=featureSubNorm2;
end

save([resultDir,'/normFeature1_', winType, '_',session,'.mat'],'normFeature1')
save([resultDir,'/normFeature2_', winType, '_',session,'.mat'],'normFeature2')

%% run PCA on standardized feature
norm='norm2';
if strcmp(norm, 'norm1')
loadFile=load([resultDir,'/normFeature1_', winType, '_',session,'.mat'])
normFeature=loadFile.normFeature1;
else 
loadFile=load([resultDir,'/normFeature2_', winType, '_',session,'.mat'])
normFeature=loadFile.normFeature2;
end

[COEFF,SCORE,latent, tsquare] = princomp(normFeature);
a=cumsum(latent)./sum(latent);
save([resultDir,'/eigenvector_',norm, '_', winType, '_',session,'.mat'],'COEFF')
save([resultDir,'/eigenvalue_',norm, '_', winType, '_',session,'.mat'],'latent')
save([resultDir,'/cumVariance_', norm, '_', winType, '_',session,'.mat'],'a')

loadFile1=load([resultDir,'/eigenvalue_norm1_', winType, '_',session,'.mat'])
eigen1=loadFile1.latent;
loadFile2=load([resultDir,'/eigenvalue_norm2_', winType, '_',session,'.mat'])
eigen2=loadFile2.latent;

loadFile3=load([resultDir,'/cumVariance_norm1_', winType, '_',session,'.mat'])
a1=loadFile3.a;
loadFile4=load([resultDir,'/cumVariance_norm2_', winType, '_',session,'.mat'])
a2=loadFile4.a;


close all
figure(1)
plot(eigen1(1:end), 'r', 'LineWidth', 2)
hold on 
plot(eigen2(1:end), 'b', 'LineWidth', 2)
xlim([-1000 16000])

figure(2)
plot(a1(1:end), 'r', 'LineWidth', 2)
hold on
plot(a2(1:end), 'b', 'LineWidth', 2)
xlim([-1000 16000])
ylim([0 1.1])
saveas(figure(1), [figDir, 'eigenValueAllWin_', session '.jpg'])
saveas(figure(2), [figDir, 'CummulVarAllWin_', session '.jpg'])


%% plot the PCs
close all

norm='norm2'
if strcmp(norm, 'norm1')
    loadFile=load([resultDir,'/eigenvectorNormFeature1_', winType, '_',session,'.mat']);
elseif strcmp(norm, 'norm2')
    loadFile=load([resultDir,'/eigenvectorNormFeature2_', winType, '_',session,'.mat']);
end
COEFF=loadFile.COEFF;

ROIOrder=xlsread('/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/doc/Cradd179_ROIorder.xlsx','Sheet1','a2:b180');
orderByLobe=ROIOrder(:,1);

close all
for i=1:12
    tmp=squareform(COEFF(:,i));
    tmpReorderCol=tmp(:, orderByLobe);
    tmpReorderColAndRow=tmpReorderCol(orderByLobe, :);
    figure(1)
    subplot(3,4,i)
    imagesc(tmpReorderColAndRow)
    colorbar
    caxis([-0.015 0.015])
    title(['PC', num2str(i)])
    
end
saveas(figure(1), [figDir, norm, 'PC_', session '.jpg'])
