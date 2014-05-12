clear
clc
close all
session='session2';
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;
numSub=length(subList)
numSeed=4;
covType='GSR'
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
corAllSub=zeros(numSeed, numSeed, numSub);
for k=1:numSub
    sub=subList{k};
    disp (['Working on sub ', char(sub),' ......'])
    subDir=[analyDir,'data/645/all_10min/',covType, '/',session,'/',char(sub)];
    seedROISignals = load([subDir,'/ROISignals_seedROISignal.mat']);
    TC=seedROISignals.ROISignals;
    numSeed1=size(TC,2);
    if numSeed1~=numSeed
        disp('numSeed dosen"t match')
    end
    corAllSub(:, :, k)=corrcoef(TC);
end

corAllSubReshape=reshape(corAllSub, [], numSub);
meanFCSeed=mean(corAllSubReshape, 2);
meanFCSeedReshape=reshape(meanFCSeed, numSeed, numSeed);

corNoDuplicate=[squeeze(corAllSub(1,2,:)) squeeze(corAllSub(1,3,:)) squeeze(corAllSub(1,4,:)) squeeze(corAllSub(2,3,:)) squeeze(corAllSub(2,4,:)) squeeze(corAllSub(3,4,:))];
[h, p]=ttest(corNoDuplicate)
[pID,pN] = FDR(p,0.05)

a=[0.200000002980232 1 1;0.257142871618271 1 1;0.314285725355148 1 1;0.371428579092026 1 1;0.428571432828903 1 1;0.485714286565781 1 1;0.54285717010498 1 1;0.600000023841858 1 1;0.657142877578735 1 1;0.714285731315613 1 1;0.77142858505249 1 1;0.828571438789368 1 1;0.885714292526245 1 1;0.942857146263123 1 1;1 1 1;1 1 0.9375;1 1 0.875;1 1 0.8125;1 1 0.75;1 1 0.6875;1 1 0.625;1 1 0.5625;1 1 0.5;1 1 0.4375;1 1 0.375;1 1 0.3125;1 1 0.25;1 1 0.1875;1 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.964705884456635 0 0;0.929411768913269 0 0;0.894117653369904 0 0;0.858823537826538 0 0;0.823529422283173 0 0;0.788235306739807 0 0;0.752941191196442 0 0;0.717647075653076 0 0;0.682352960109711 0 0;0.647058844566345 0 0;0.61176472902298 0 0;0.576470613479614 0 0;0.541176497936249 0 0;0.505882382392883 0 0;0.470588237047195 0 0;0.43529412150383 0 0;0.400000005960464 0 0];
figure(1)
imagesc(meanFCSeedReshape)
caxis([-0.1 0.4])

set(gca, 'YTick', [], 'XTick', [])
% title('Mean')
% xlabel('Seeds')
% ylabel('Seeds')
%colormap(cmap)
colormap(a)
hcb=colorbar
set(hcb,'YTick',[-0.1, 0, 0.1, 0.2, 0.3, 0.4])
saveas(figure(1), sprintf('/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/GSR/%s/stationaryFCSeed_%s.png', session, session))
