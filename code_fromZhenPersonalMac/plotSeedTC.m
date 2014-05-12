%This script compute thess average sateionary FC across all
% subjects
clear
clc
close all


sessionList={'session1','session2'};
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;
numSeed=4;
numROI=156;


%sessionList={'session1'};
% subList={'0021002'};

numSub=length(subList);
numSession=length(sessionList);
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
figDir=[analyDir, 'fig/645/all_10min/seedTC/']

% The time courses of four seeds were plotted for ession 1, sub 3 and 14 and session 2, sub21
for j=1
    session=char(sessionList{j});
    resultDir=[analyDir,'results/645/all_10min/',session,'/'];
    for k=21
        sub=subList{k};
        disp (['Working on sub ', char(sub),' ......'])
        subDir=[analyDir,'data/645/all_10min/',session,'/',char(sub)];
        seedROISignals = load([subDir,'/ROISignals_seedROISignal.mat']);
        TC1=seedROISignals.ROISignals;
        figure(1)
        for i=1:numSeed
            subplot(2,2,i)
            plot(TC1(:,i))
            xlabel('Time')
            ylabel('Intensity')
        end
%         ps=ginput(1)
%         text(ps(1),ps(2),[session, 'sub',sub, ' Seed Time Course'])
        saveas(figure(1),[figDir, 'seedTC_',session,'_',sub,'.png'])
    end
end



