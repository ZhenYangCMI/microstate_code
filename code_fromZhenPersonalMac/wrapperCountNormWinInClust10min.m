
clear
clc
close all

sessionList={'session1','session2'};
%sessionList={'session1'}
numSession=length(sessionList);
covTypeList={'Full'};
%covTypeList={'Full','Partial'};
numCovType=length(covTypeList);
numSeed=4;
dataLength='all_10min';

q=0;
p=0;
for i=1:numSession
    session=char(sessionList{i});
    dataDir= (['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength, filesep, session,'/']);
    figDir=(['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength,filesep, session,'/']);
    
    for j=1:numCovType
        covType=char(covTypeList{j});
        q=q+1;
        % coumpute and plot number of windows in a cluster for clustering
        % with all seeds concatenated
        indxAllSeeds=load([dataDir,'zWinAllSeeds_',covType, 'CorLasso_',session, '_10min.txt']);
        [ numWinInClustAllSeeds ] = countWinInClust(indxAllSeeds);
        numClust1(q)=length(numWinInClustAllSeeds);
        clustNumAllSeeds=1:length(numWinInClustAllSeeds);
        figure(q)
        bar(clustNumAllSeeds,numWinInClustAllSeeds)
        title(['Number of Nomralized Windows in a Cluster All Seeds ',covType, 'Cor ', session])
        xlabel('Cluster Number');
        ylabel('Number of Windows');
        if strcmp(covType,'Full')
            set(gca,'XTick',0:7)
            xlim([0 7])
            set(gca,'YTick',0:1000:6000);
            ylim([0 6000]);
        else
            set(gca,'XTick',0:25)
            xlim([0 25])
            set(gca,'YTick',0:500:4000);
            ylim([0 4000]);
        end
        saveas(figure(q),[figDir,'zWinAllSeeds_',covType, 'CorLasso_',session,'.png'])
        
        
        % coumpute and plot number of windows in a cluster for clustering
        % with windows of one seed
        indx=load([dataDir,'zWinEachSeed_',covType, 'CorLasso_',session, '_10min.txt']);
        for k=1:numSeed
            p=p+1;
            indxEachSeed=squeeze(indx(:,k));
            [ numWinInClustEachSeed] = countWinInClust(indxEachSeed);
            numClust2(p)=length(numWinInClustEachSeed);
            clustNumEachSeed=1:length(numWinInClustEachSeed);
            figure(q+4)
            subplot(2,2,k)
            bar(clustNumEachSeed,numWinInClustEachSeed)
            title(['seed ',num2str(k)])
            xlabel('Cluster Number');
            ylabel('Number of Windows');
            if strcmp(covType,'Full')
                set(gca,'XTick',0:2:14)
                xlim([0 14])
                set(gca,'YTick',0:300:1800);
                ylim([0 1800]);
            else
                set(gca,'XTick',0:15:60)
                xlim([0 60])
                set(gca,'YTick',0:300:1500);
                ylim([0 1500]);
            end
            hold on
        end
        hold off
        saveas(figure(q+4),[figDir,'zWinEachSeed_',covType, 'CorLasso_',session,'.png'])
    end
end