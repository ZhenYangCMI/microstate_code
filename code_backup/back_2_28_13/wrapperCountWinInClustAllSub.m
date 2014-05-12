
clear
clc
close all

%sessionList={'session1','session2'};
sessionList={'session1'}
numSession=length(sessionList);
covTypeList={'Partial'};
%covTypeList={'Full','Partial'};
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
        indxAllSeeds=load([dataDir,'winAllSeeds_',covType, 'CorLasso_',session, '.txt']);
        [ numWinInClustAllSeeds ] = countWinInClust(indxAllSeeds);
        numClust1(q)=length(numWinInClustAllSeeds);
        clustNumAllSeeds=1:length(numWinInClustAllSeeds);
        figure(q)
        bar(clustNumAllSeeds,numWinInClustAllSeeds)
        title(['Number of Windows in a Cluster All Seeds ',covType, 'Cor ', session])
        xlabel('Cluster Number');
        ylabel('Number of Windows');
        if strcmp(covType,'Full')
            set(gca,'XTick',0:9)
            xlim([0 9])
            set(gca,'YTick',0:500:4000);
            ylim([0 4000]);
        else
            set(gca,'XTick',0:25)
            xlim([0 25])
            set(gca,'YTick',0:500:4000);
            ylim([0 4000]);
        end
        saveas(figure(q),[figDir,'winAllSeeds_',covType, 'CorLasso_',session,'.png'])
        
        
        % coumpute and plot number of windows in a cluster for clustering
        % with windows of one seed
        indx=load([dataDir,'winEachSeed_',covType, 'CorLasso_',session, '.txt']);
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
                set(gca,'XTick',0:14)
                xlim([0 14])
                set(gca,'YTick',0:300:1500);
                ylim([0 1500]);
            else
                set(gca,'XTick',0:15:60)
                xlim([0 60])
                set(gca,'YTick',0:300:1500);
                ylim([0 1500]);
            end
            hold on
        end
        hold off
        saveas(figure(q+4),[figDir,'winEachSeed_',covType, 'CorLasso_',session,'.png'])
    end
end