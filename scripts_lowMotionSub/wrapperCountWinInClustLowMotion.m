
clear
clc
close all

sessionList={'session1','session2'};
numSession=length(sessionList);
%covTypeList={'Full'};
covTypeList={'Full','Partial'};
numCovType=length(covTypeList);
numSeed=4;

q=0;
p=0;
for i=1:numSession
    session=char(sessionList{i});
    dataDir= (['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',session,'/',session,'_ward_euclidean/lambdaOptimalPerSub/']);
    figDir=(['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',session,'/',session,'_ward_euclidean/lambdaOptimalPerSub/']);
    
    for j=1:numCovType
        covType=char(covTypeList{j});
        q=q+1;
        % coumpute and plot number of windows in a cluster for clustering
        % with all seeds concatenated
        indxAllSeeds=load([dataDir,'winAllSeeds_',covType, 'CorLasso_',session, '_lowMot.txt']);
        [ numWinInClustAllSeeds ] = countWinInClust(indxAllSeeds);
        numClust1(q)=length(numWinInClustAllSeeds);
        clustNumAllSeeds=1:length(numWinInClustAllSeeds);
        figure(q)
        bar(clustNumAllSeeds,numWinInClustAllSeeds)
        title(['Number of Windows in a Cluster All Seeds ',covType, 'Cor ', session])
        xlabel('Cluster Number');
        ylabel('Number of Windows');
        if strcmp(covType,'Full')
            set(gca,'XTick',0:7)
            xlim([0 7])
            set(gca,'YTick',0:200:1200);
            ylim([0 1200]);
        else
            set(gca,'XTick',0:5:45)
            xlim([0 45])
            set(gca,'YTick',0:100:600);
            ylim([0 600]);
        end
        saveas(figure(q),[figDir,'winAllSeeds_',covType, 'CorLasso_',session,'_LowMotion.png'])
        
        
        % coumpute and plot number of windows in a cluster for clustering
        % with windows of one seed
        indx=load([dataDir,'winEachSeed_',covType, 'CorLasso_',session, '_lowMot.txt']);
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
                set(gca,'XTick',0:2:20)
                xlim([0 20])
                set(gca,'YTick',0:50:300);
                ylim([0 300]);
            else
                set(gca,'XTick',0:10:40)
                xlim([0 40])
                set(gca,'YTick',0:30:150);
                ylim([0 150]);
            end
            hold on
        end
        hold off
        saveas(figure(q+4),[figDir,'winEachSeed_',covType, 'CorLasso_',session,'_LowMotion.png'])
    end
end