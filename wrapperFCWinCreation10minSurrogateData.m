%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script generate full and partial correlation windows for all
% subjects
clear
clc
close all


parfor surrogateNum=51:100 % num of surrogate
    disp(['Working on surrogate_', num2str(surrogateNum)])
    %TRList={'645','2500'};
    %sessionList={'session1','session2'};
    subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
    subList=subList.SubID;
    
    TRList={'645'};
    sessionList={'session1'};
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
            resultDir=[analyDir,'results/',TR, '/all_10min/',session,'/surrogate/'];
            tmp=load([analyDir,'results/645/all_10min/meanLambdaAtMaxLogAcrossWin.txt']);
            lambdaList=squeeze(tmp(:,j));
            %clear winFullCor winFullCorLasso winPartialCorLasso
            [winFullCor, winFullCorLasso, winPartialCorLasso]=FCWinCreationSurrogateDataAllSub(TR,session,subList, lambdaList, surrogateNum);
            %save([resultDir,'winFullCor_Corrcoef_',TR,'_',session,'.mat'],'winFullCor')
            %save([resultDir,'winFullCorLasso_OptimalLambdaPerSub_', TR,'_',session,'_surrogate_', num2str(k),'.mat'],'winFullCorLasso')
            %save([resultDir,'winPartialCorLasso_OptimalLambdaPerSub_', TR,'_',session,'.mat'],'winPartialCorLasso')
        end
    end
    
end