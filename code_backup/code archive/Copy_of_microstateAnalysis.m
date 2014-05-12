% This script will run the microstate analysis

clear
clc
close all

%TRList={'645','2500'};
sessionList={'session1','session2'};
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;

TRList={'645'};
%sessionList={'session1'};
%subList={'0021002'};

numSub=length(subList);
numSession=length(sessionList);
numTR=length(TRList);
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];


% generate full and partial correlation windows
for i=1:numTR
    TR=TRList{i};
    for j=1:numSession
        session=sessionList{j};
        resultDir=[analyDir,'/results/',TR, '/',session];
        clear winAllSubAllSeed partWinAllSubAllSeed
        [winAllSubAllSeedFullCor, winAllSubAllSeedPartialCor]=FCWinCreation(TR,session,subList);
        save([resultDir,'/winAllSubAllSeed_fullCor_',char(TR),'_',char(session),'.mat'],'winAllSubAllSeedFullCor');
        save([resultDir,'/winAllSubAllSeed_partialCor_',char(TR),'_',char(session),'.mat'],'WinAllSubAllSeedPartialCor');
    end
end

% estimate cluster numbers
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
        
        corType=['full'];
        cutoffThresh=0.1;
        clear indx numclust clustTree cophenetCoef inconsistCoef
        [indx,numclust,clustTree,cophenetCoef,inconsistCoef]=hierarchicalCluster...
            (corType, cutoffThresh, winAllSubAllSeedFullCor, winAllSubAllSeedPartialCor);
        save([resultDir,'/fullCorIndx_',TR,'_',session],'indx')
        save([resultDir,'/fullCorAllVariables_',TR,'_',session])
        figure(1)
        dendrogram(clustTree)
        saveas(figure(1),[figDir,'/dendrogram_',TR,'_',session,'_',corType,'.png'])
        disp('Hierarchical clustering of the full correlation windows done.')
        
        corType=['partial'];
        cutoffThresh=0.1;
        clear indx numclust clustTree cophenetCoef inconsistCoef
        [indx,numclust,clustTree,cophenetCoef,inconsistCoef]=hierarchicalCluster...
            (corType, cutoffThresh, winAllSubAllSeedFullCor, winAllSubAllSeedPartialCor);
        save([resultDir,'/partialCorIndx_',TR,'_',session],'indx')
        save([resultDir,'/partialCorAllVariables_',TR,'_',session])
        figure(2)
        dendrogram(clustTree)
        saveas(figure(2),[figDir,'/dendrogram_',TR,'_',session,'_',corType,'.png'])
        disp('Hierarchical clustering of the partial correlation windows done.')
        
    end
end


