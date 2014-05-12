clear
clc
close all

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
        
corWin=corrcoef( winAllSubAllSeedFullCor');
                imagesc(corWin)