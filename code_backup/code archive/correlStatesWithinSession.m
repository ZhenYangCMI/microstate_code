clear
clc

resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/session2/ses2_ward_euclidean/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/session2/ses2_ward_euclidean/'];

k=0;

for i=2:2:16
    k=k+1;
    numClust=i;
    tmp=load([resultDir,'clusterMean_',num2str(numClust),'clusters.mat']);
    clustMean=tmp.finalMeanWinOfClust;
    clustMeanTransp=clustMean';
    corClusters=corrcoef(clustMeanTransp); 
        
    figure(k)
    imagesc(corClusters)
    colorbar
    title('Correlations between clusters')
    xlabel('Clusters')
    ylabel('clusters')
    set(gca,'XTick',-1:1,'YTick',-1:1);
    saveas(figure(k),[figDir,'Correlations_',num2str(numClust),'clusters.png'])
end