% This script generate full and partial correlation windows
clear
clc
close all

%TRList={'645','2500'};
%sessionList={'session1','session2'};
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;

TRList={'645'};
sessionList={'session2'};
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
        save([resultDir,'/winAllSubAllSeed_fullCor_',char(TR),'_',char(session),'_0.12.mat'],'winAllSubAllSeedFullCor');
        save([resultDir,'/winAllSubAllSeed_partialCor_',char(TR),'_',char(session),'_0.12.mat'],'winAllSubAllSeedPartialCor');
    end
end