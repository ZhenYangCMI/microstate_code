clear
clc
close all

% Load the full and partial correlation windows
TRList={'645'};
sessionList={'session2'};
numSession=length(sessionList);
numTR=length(TRList);
linkageMethodList={'ward','average','centroid','complete','median','single','weighted'};
distMetricList={'euclidean','seuclidean','cityblock','minkowski','chebychev','mahalanobis','cosine','correlation','spearman','hamming','jaccard'};


numLinkageMethod=length(linkageMethodList);
numDistMetric=length(distMetricList);

analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed'];

for i=1:numTR
    TR=TRList{i};
    for j=1:numSession
        session=sessionList{j};
        disp(['Working on ',TR,' ', session, '.'])
        resultDir=[analyDir,'/results/',TR, '/',session];
        figDir=[analyDir,'/fig/',TR,'/',session];
        tmp1=load([resultDir,'/winAllSubAllSeed_fullCor_',char(TR),'_',char(session),'.mat']);
        tmp2=load([resultDir,'/winAllSubAllSeed_partialCor_',char(TR),'_',char(session),'.mat']);
        winAllSubAllSeedFullCor=tmp1.winAllSubAllSeed;
        winAllSubAllSeedPartialCor=tmp2.partWinAllSubAllSeed;
        disp('Windows were loaded successfully.')
    end
end

% estimate different linkage methods and distance metrics


corType=['partial'];
if strcmp('full',corType)
    winForClustering=winAllSubAllSeedFullCor;
elseif strcmp('partial',corType)
    winForClustering=winAllSubAllSeedPartialCor;
else
    disp('Incorrect input for corType.')
end

inconsistCoef=zeros(numWin,4,numLindageMethod*numDistMetric);
k=0;
for m=6:numLinkageMethod
    linkageMethod=linkageMethodList{m};
    for n=1:numDistMetric
        distMetric=distMetricList{n};
        k=k+1
        clear distBtwWin clustTree inconsistCoef
        disp(['Working on ',linkageMethod,'_',distMetric,'.'])
        distBtwWin= pdist(winForClustering,distMetric);
        clustTree = linkage(distBtwWin,linkageMethod);
        cophenetCoef(k) = cophenet(clustTree,distBtwWin)
        inconsistCoef=inconsistent(clustTree);
        figure(k)
        dendrogram(clustTree,0)
        title(['Dentrogram ',linkageMethod,'',distMetric])
        saveas(figure(k),[figDir,'/dendrogram_',corType,'_',linkageMethod,'_',distMetric,'.png'])
    end
end
save([resultDir,'/',TR,'_',session,'_',corType,'_cophenetCoef.mat'],'cophenetCoef')




