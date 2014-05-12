clear
clc
close all


dataLength='all_10min';

resultDir1=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/', dataLength, filesep, 'GSR/session1/clustMean/'];
resultDir2=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/', dataLength, filesep, 'GSR/session2/clustMean/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/', dataLength, '/correlTwoSessions/'];


    
    numClust1=5;
    numClust2=6;
    tmp1=load([resultDir1,'/clusterMean_',num2str(numClust1),'clusters_session1_normWin.mat']);
    tmp2=load([resultDir2,'/clusterMean_',num2str(numClust2),'clusters_session2_normWin.mat']);
    clustMean1=tmp1.finalMeanWinOfClust;
    tmp3=clustMean1';
    clustMean2=tmp2.finalMeanWinOfClust;
    tmp4=clustMean2';
    
    % plot the distribution of the clustMeans
    figure(1)
for i=1:5
    if i<=5
        subplot(2,3,i)
        hist(clustMean1(i,:))
        ylim([0 40])
    end
end

figure(2)
for i=1:6
    subplot(2,3,i)
    hist(clustMean2(i,:))
    ylim([0 40])
end

    
    
    % recoder the cluster in session 2 to make them match with session
    % according to correlation value
    clustMeanTransp1=horzcat(tmp3(:,2), tmp3(:,1), tmp3(:,3), tmp3(:, 4),tmp3(:,5));
    clustMeanTransp2=horzcat(tmp4(:,2), tmp4(:,6), tmp4(:,4), tmp4(:, 1), tmp4(:,5), tmp4(:,3));
    
    clustMeanTransp=[clustMeanTransp1,clustMeanTransp2];
    [corClusters,pValue]=corrcoef(clustMeanTransp);
    
    [pID,pN] = FDR(pValue(6:end, 1:5),0.05)
    pValue(6:end, 1:5)
    
    corPlot=corClusters(6:end, 1:5);
    load('MyColormapsCorrel','mycmap')
    
    figure(1)
    imagesc(corPlot)
    h=figure(1);
    %colorbar
    caxis([-1 1])
    set(gca,'xTick',0.5:5.5, 'yTick', 0.5:5.5,'XTickLabel',[],'YTickLabel',[]);
    set(gca,'LineWidth',2)
    set(figure(1),'Colormap',mycmap)
    %set(figure(1),'Colormap',ColorMap)  % this is used to generate the
    %first graph for eiditing the colormap
    grid on
    set(gca, 'GridLineStyle', '-' )
    
    set(h,'Position',[1,1,800,665])
    saveas(figure(1),[figDir,'CorrelSes1Ses2_zWinFullCorLasso.png'])
    
    
    figure(1)
    imagesc(corClusters)
    colorbar
    title('Correlations between clusters')
    xlabel('Clusters')
    ylabel('clusters')
    caxis([-1 1])
    saveas(figure(1),[figDir,'CorrelBetwTwoSessions_zWinFullCorLasso.png'])



