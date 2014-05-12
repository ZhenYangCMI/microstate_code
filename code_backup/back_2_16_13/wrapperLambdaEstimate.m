clear
clc
close all

%TRList={'645','2500'};
sessionList={'session1','session2'};
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;

TRList={'645'};
% sessionList={'session1'};
% subList={'0021002','0021018'};
winWidth=69;
step=69;
TR=0.645;
numVol=450;
numWin=floor((numVol-winWidth)/step)+1;

numSub=length(subList);
numSession=length(sessionList);
numTR=length(TRList);
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
figDir=[analyDir,'fig/lambdaEstimate/'];
resultDir=[analyDir,'results/645'];

lambdaAtMaxLog=zeros(numSub,numWin,numSession);
for i=1:numSession
    session=char(sessionList{i});
    
    for j=1:numSub
        subID=subList{j};
        disp(['Working on ',session, '_',subID, '.'])
        [logLikelihoodList, lambdaList]=lambdaEstimate(session, subID);
        disp(['LogLikelihoodList and lambdaList for ', session, '_',subID, 'computed.'])
        for k=1:numWin
            logLikelihood=logLikelihoodList(k,:);
            lambda=lambdaList(k,:);
            Max=max(logLikelihood);
            lambdaAtMaxLog(j,k,i)=lambda(find(logLikelihood==Max));
        end
    end
end
save([resultDir,'/lambdaAtMaxLog_',num2str(numSession),'sessions_',num2str(numSub),'subjects.mat'], 'lambdaAtMaxLog')
disp('lambdaAtMaxLog fro all subjects and all sessions found and saved in one matrix.')


% plot the correlations between windows and subjects
lambdaAtMaxLog2D=reshape(lambdaAtMaxLog,numSub,[]);

corWinAllSession=corrcoef(lambdaAtMaxLog2D);
figure(1)
imagesc(corWinAllSession)
title('corWinAllSession')
colorbar('EastOutside')
ylabel('LambdaAtMaxLog for Each Window')
xlabel('LambdaAtMaxLog for Each Window')
saveas(figure(1),[figDir, 'corWinAllSession.png']);

corSub=corrcoef(lambdaAtMaxLog2D');
figure(2)
imagesc(corSub)
title('corSub')
colorbar('EastOutside')
ylabel('LambdaAtMaxLog for Each Subject')
xlabel('LambdaAtMaxLog for Each Subject')
saveas(figure(2),[figDir, 'corSub.png']);


% plot the distribution of lambda
lambdaAtMaxLog1D=reshape(lambdaAtMaxLog,1,[]);
figure(3)
subplot(1,2,1)
plot(lambdaAtMaxLog1D)
ylabel('Lambda Value')
xlabel('Number of Lambda')
subplot(1,2,2)
hist(lambdaAtMaxLog1D,0.05:0.01:0.2)
ylabel('Frequency')
xlabel('Lambda Value')
saveas(figure(3),[figDir, 'Lambda Distribution.png']);

