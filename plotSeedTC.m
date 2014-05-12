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

% plot time course for presentation
sub2Plot=[10];

p=0;
for j=1:numSession
    session=char(sessionList{j});
    resultDir=[analyDir,'results/645/all_10min/',session,'/'];
    for k=1:length(sub2Plot)
        indx=sub2Plot(k);
        sub=subList{indx};
        disp (['Working on sub ', char(sub),' ......'])
        subDir=[analyDir,'data/645/all_10min/',session,'/',char(sub)];
        seedROISignals = load([subDir,'/ROISignals_seedROISignal.mat']);
        TC1=seedROISignals.ROISignals;
        p=p+1;
        figure(p)
        for i=1:numSeed
            x=TC1(:,i);
            y=x(1:3:end);
            plot(y, 'Color',[0.8-0.2*i 0.8-0.2*i 0.8-0.2*i],'LineWidth',0.2)
                                  ylim([-5.5 5.5])
                                  axis on
            hold on
        end
        %saveas(figure(p),[figDir, 'seedTC_',session,'_',sub,'.png'])
    end
end



% The time courses of four seeds were plotted for ession 1, sub 3 and 14 and session 2, sub21
sub2Plot=[3,14,21];
p=0;
for j=1:numSession
    session=char(sessionList{j});
    resultDir=[analyDir,'results/645/all_10min/',session,'/'];
    for k=1:length(sub2Plot)
        indx=sub2Plot(k);
        sub=subList{indx};
        disp (['Working on sub ', char(sub),' ......'])
        subDir=[analyDir,'data/645/all_10min/',session,'/',char(sub)];
        seedROISignals = load([subDir,'/ROISignals_seedROISignal.mat']);
        TC1=seedROISignals.ROISignals;
        p=p+1;
        figure(p)
        for i=1:numSeed
            subplot(2,2,i)
            plot(TC1(:,i))
            xlabel('Time')
            ylabel('Intensity')
            set(gca, 'YTick', [-4 -2 0 2 4])
            ylim([-4.5 4.5])
            title(['Seed ', num2str(i)])
        end
        text=['Seed Time Courses: sub', sub, session];
        ht = subtitle(2,text)
        saveas(figure(p),[figDir, 'seedTC_',session,'_',sub,'.png'])
    end
end



