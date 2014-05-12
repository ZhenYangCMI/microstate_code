%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script will plot the distribution of all features for each type of window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

% load subList
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;
%subList={'0021002', '0021006'};
numSub=length(subList);

session='session1'

% define the analysis strategy, use to fid the right TS data
covType='GSR';

% define the windodw parameters
winSize=69;
step=3; % in TR
numVol=884;

%define data and resultDir
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
dataDir=[analyDir, '/data/645/all_10min/', covType, '/',session,'/'];
resultDir=[analyDir,'fullBrainFC/results/',covType, '/', session,'/'];
figDir=[analyDir, 'fullBrainFC/figs/GSR/cmpWinType/']

% 1 Extract the features from the full brain correlation matrix and
% standardize the features for each featureWin

tmp1=load([resultDir,'featureFC_winFullCor_', session,'.mat'])
tmp2=load([resultDir,'featureFC_winFullCorLasso_',session,'.mat'])
tmp3=load([resultDir,'featureFC_winPartialCorLasso_', session,'.mat'])
winFullCor=tmp1.featureWin;
winFullCorLasso=tmp2.featureWin;
winPartialCorLasso=tmp3.featureWin;

% plot the correlation distribution
r1=(exp(2*winFullCor)-1)./(exp(2*winFullCor)+1);
r2=(exp(2*winFullCorLasso)-1)./(exp(2*winFullCorLasso)+1);
r3=-(exp(2*winPartialCorLasso)-1)./(exp(2*winPartialCorLasso)+1);

r1Reshap=reshape(r1, [],1);
r2Reshap=reshape(r2, [],1);
r3Reshap=reshape(r3, [],1);
r=horzcat(r1Reshap, r2Reshap, r3Reshap);
winType={'winFullCor', 'winFullCorLasso', 'winPartialCorLasso'};
figure(1)
for j=1:length(winType)
        subplot(1,3,j)
    win=char(winType{j})
    xvalue=-1.5:0.25:1.5;
    
    if j==1
        hist(r(:,1), xvalue)
    else
        clear h
        hist(r(:,1), xvalue)
        hold on
        hist(r(:,j), xvalue)
        h = findobj(gca,'Type','patch');
        set(h(1),'FaceColor','r', 'EdgeColor','k')
        set (h(1), 'FaceAlpha', 0.3);
    end
    title(win)
    ylim([0 100000000])
end
saveas(figure(1), [figDir, session, '_cmpWinTypeCorDistr.png'])








