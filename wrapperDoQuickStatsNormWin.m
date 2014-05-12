
clear
clc
close all

numSeed=4;
numSub=22;
sessionList={'session1','session2'};
numSession=length(sessionList);

for i=1:numSession
    close all
    session=char(sessionList{i});
    resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',session,'/',session,'_ward_euclidean/lambdaOptimalPerSub/'];
    figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',session,'/',session,'_ward_euclidean/lambdaOptimalPerSub/percent/'];
    
    
    indx=load([resultDir,'zWinAllSeeds_FullCorLasso_',session,'.txt']);
    numClust=length(unique(indx));
    [seedPrct, subPrct, subSeedPrct] = DoQuickStats(indx,numClust);
    
    % plot seedPrct
    figure(9)
    imagesc(seedPrct')
    colorbar
    title('Percent of Norm windows in each seed across sub')
    ylabel('Seeds')
    xlabel('clusters')
    caxis([0 1])
    %     set(gca,'XTick',0:9)
    %             xlim([0 9])
    set(gca,'YTick',1:4);
    ylim([0.5 4.5]);
    saveas(figure(9),[figDir,'seedPrct_',num2str(numClust),'clusters_normWin.png'])
    
    % plot subPrct
    figure(10)
    imagesc(subPrct')
    colorbar
    title('Percent of Norm windows in each sub across seed')
    ylabel('Subjects')
    xlabel('clusters')
    caxis([0 1])
    saveas(figure(10),[figDir,'subPrct_',num2str(numClust),'clusters_normWin.png'])
    
    % plot subSeedPrct
    categoryName='subSeedPrct';
    for m=1:numSub
        clear dataPlot1
        dataPlot1=squeeze(subSeedPrct(:,:,m))';
        if (m<5)
            figure(1)
            subplot(2,2,m)
            imagesc(dataPlot1)
            title('Sub 1 to sub 4')
            colorbar
            ylabel('Seeds')
            xlabel('clusters')
            caxis([0 1])
        elseif (m<9)
            figure(2)
            subplot(2,2,m-4)
            imagesc(dataPlot1)
            title('Sub 5 to sub 8')
            colorbar
            ylabel('Seeds')
            xlabel('clusters')
            caxis([0 1])
        elseif (m<13)
            figure(3)
            subplot(2,2,m-8)
            imagesc(dataPlot1)
            title('Sub 9 to sub 12')
            colorbar
            ylabel('Seeds')
            xlabel('clusters')
            caxis([0 1])
        elseif (m<17)
            figure(4)
            subplot(2,2,m-12)
            imagesc(dataPlot1)
            title('Sub 13 to sub 16')
            colorbar
            ylabel('Seeds')
            xlabel('clusters')
            caxis([0 1])
        elseif (m<21)
            figure(5)
            subplot(2,2,m-16)
            imagesc(dataPlot1)
            title('Sub 17 to sub 20')
            colorbar
            ylabel('Seeds')
            xlabel('clusters')
            caxis([0 1])
        else
            figure(6)
            subplot(2,2,m-20)
            imagesc(dataPlot1)
            title('Sub 21 to sub 22')
            colorbar
            ylabel('Seeds')
            xlabel('clusters')
            caxis([0 1])
        end
    end
    saveas(figure(1),[figDir,categoryName,'_plotPerSeed_1-4Sub_',num2str(numClust),'clusters_normWin.png'])
    saveas(figure(2),[figDir,categoryName,'_plotPerSeed_5-8Sub_',num2str(numClust),'clusters_normWin.png'])
    saveas(figure(3),[figDir,categoryName,'_plotPerSeed_9-12sub_',num2str(numClust),'clusters_normWin.png'])
    saveas(figure(4),[figDir,categoryName,'_plotPerSeed_13-16Sub_',num2str(numClust),'clusters_normWin.png'])
    saveas(figure(5),[figDir,categoryName,'_plotPerSeed_17-20Sub_',num2str(numClust),'clusters_normWin.png'])
    saveas(figure(6),[figDir,categoryName,'_plotPerSeed_21-22sub_',num2str(numClust),'clusters_normWin.png'])
    
    for n=1:numSeed
        clear dataPlot2
        dataPlot2=squeeze(subSeedPrct(:,n,:))';
        
        if (n<3)
            figure(7)
            subplot(1,2,n)
            imagesc(dataPlot2)
            colorbar
            title('Seed 1 and Seed 2')
            ylabel('Subjects')
            xlabel('clusters')
            caxis([0 1])
        else
            figure(8)
            subplot(1,2,n-2)
            imagesc(dataPlot2)
            colorbar
            title('Seed 3 and Seed 4')
            ylabel('Subjects')
            xlabel('clusters')
            caxis([0 1])
        end
        saveas(figure(7),[figDir,categoryName,'_plotPerSub_first2Seeds_',num2str(numClust),'clusters_normWin.png'])
        saveas(figure(8),[figDir,categoryName,'_plotPerSub_last2Seeds_',num2str(numClust),'clusters_normWin.png'])
    end
    
end

