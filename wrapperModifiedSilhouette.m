clear
clc
close all

%TRList={'645','2500'};
%sessionList={'session1','session2'};

TRList={'645'};
sessionList={'session1'};

numSession=length(sessionList);
numTR=length(TRList);
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];


minNumClust=2;
maxNumClust=100;
numEstimate=maxNumClust-minNumClust+1;
for i=1:numTR
    TR=TRList{i};
    for j=1:numSession
        session=sessionList{j};
        resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min',filesep, session,'/'];
        figDir=[analyDir,'/fig/',TR,'/all_10min', filesep, session,'/'];
        tmp1=load([resultDir,'zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'.mat']);
        zWinFullCorLasso=tmp1.zWinFullCorLasso;
        similarityMax=corrcoef(zWinFullCorLasso');
        indxList=load([resultDir,'clustIndx_',num2str(minNumClust),'to',num2str(maxNumClust), 'clusters_zWinAllSeeds_FullCorLasso_', session, '_10min.txt']);
        
        disp ('The window and index files are loaded successfully.')
        
        msilhouette=zeros(numEstimate)
        for k=minNumClust:maxNumClust
            [msilhouette(k-1)] = modified_silhouette(similarityMat, k, squeeze(indxList(:,k-1)))
        end
        figure (1)
        plot(msilhouette)
        title('Modified Silhouette as a function of number of clusters')
        xlabel('Number of Clusters')
        ylabel('Modified Silhouette')
        saveas(figure(1),[figDir,'mSilhouette_',session,'_',num2str(minNumClust),'_to_', num2str(maxNumClust),'clusters.png'])
    end
end