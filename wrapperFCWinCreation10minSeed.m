%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script generate full and partial correlation windows for all
% subjects
clear
clc
close all

%TRList={'645','2500'};
%sessionList={'session1','session2'};
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;

TRList={'645'};
sessionList={'session2'};
% subList={'0021002'};
covType='GSR';
winSize=69;
numSub=length(subList);
numSession=length(sessionList);
numTR=length(TRList);
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];


% generate full and partial correlation windows
for i=1:numTR
    TR=char(TRList{i});
    for j=1:numSession
        session=char(sessionList{j});
        resultDir=[analyDir,'results/',TR, '/all_10min/',covType, filesep, session,'/'];
        tmp=load([analyDir,'results/645/all_10min/meanLambdaAtMaxLogAcrossWin.txt']);
        if strcmp(session, 'session1')
        lambdaList=squeeze(tmp(:,1));
        else
            lambdaList=squeeze(tmp(:,2));
        end
        clear winFullCorSeed winFullCorLassoSeed winPartialCorLassoSeed
        [winFullCorSeed, winFullCorLassoSeed, winPartialCorLassoSeed]=FCWinCreationBtwSeeds(TR,session,subList, lambdaList, covType);
        save([resultDir,'winFullCorSeed_Corrcoef_',TR,'_',session,'.mat'],'winFullCorSeed')
        save([resultDir,'winFullCorLassoSeed_OptimalLambdaPerSub_', TR,'_',session,'win', num2str(winSize), '_',covType,'.mat'],'winFullCorLassoSeed')
        save([resultDir,'winPartialCorLassoSeed_OptimalLambdaPerSub_', TR,'_',session,'.mat'],'winPartialCorLassoSeed')
    end
end
