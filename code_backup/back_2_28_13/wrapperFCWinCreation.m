%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script generate full and partial correlation windows for all
% subjects
clear
clc
close all

%TRList={'645','2500'};
sessionList={'session1','session2'};
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;

TRList={'645'};
% sessionList={'session2'};
% subList={'0021002'};

numSub=length(subList);
numSession=length(sessionList);
numTR=length(TRList);
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];


% generate full and partial correlation windows
for i=1:numTR
    TR=char(TRList{i});
    for j=1:numSession
        session=char(sessionList{j});
        resultDir=[analyDir,'results/',TR, '/',session,'/'];
        tmp=load([analyDir,'results/645/meanLambdaAtMaxLogAcrossWin.txt']);
        lambdaList=squeeze(tmp(:,j));
        clear winFullCor winFullCorLasso winPartialCorLasso
        [winFullCor, winFullCorLasso, winPartialCorLasso]=FCWinCreation(TR,session,subList, lambdaList);
        save([resultDir,'winFullCor_Corrcoef_',TR,'_',session,'.mat'],'winFullCor')
        save([resultDir,'winFullCorLasso_OptimalLambdaPerSub_', TR,'_',session,'.mat'],'winFullCorLasso')
        save([resultDir,'winPartialCorLasso_OptimalLambdaPerSub_', TR,'_',session,'.mat'],'winPartialCorLasso')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script generate full and partial correlation windows for one tercile
%of the subjects with low motion
clear
clc
close all

%TRList={'645','2500'};
sessionList={'session1','session2'};

TRList={'645'};
% sessionList={'session2'};
% subList={'0021002'};

numSession=length(sessionList);
numTR=length(TRList);
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];


% generate full and partial correlation windows
for i=1:numTR
    TR=char(TRList{i});
    for j=1:numSession
        session=char(sessionList{j});
        subList=textread(['/home/data/Projects/microstate/subListLowMotion',session,'.txt'],'%s');
        resultDir=[analyDir,'results/',TR, '/',session,'/'];
        lambdaList=load(['/home/data/Projects/microstate/lambdaAtMaxLogLowMotion',session,'.txt']);
        clear winFullCor winFullCorLasso winPartialCorLasso
        [winFullCor, winFullCorLasso, winPartialCorLasso]=FCWinCreation(TR,session,subList, lambdaList);
        save([resultDir,'winFullCor_Corrcoef_',TR,'_',session,'_LowMotion.mat'],'winFullCor')
        save([resultDir,'winFullCorLasso_OptimalLambdaPerSub_', TR,'_',session,'_LowMotion.mat'],'winFullCorLasso')
        save([resultDir,'winPartialCorLasso_OptimalLambdaPerSub_', TR,'_',session,'_LowMotion.mat'],'winPartialCorLasso')
    end
end


