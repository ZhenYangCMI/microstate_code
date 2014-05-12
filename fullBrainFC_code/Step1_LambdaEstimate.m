%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script estimate the optimum lambda for each subject based on the maximum log likelihood
%If two sessions were run, the lambdas for two sessions will be saved in one .mat file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

covType='GSR';

% Load the subjectList
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;
% subList={'0021002','0021018'};


% define the session
%sessionList={'session1','session2'};
sessionList={'session2'};

% define the window parameters
winSize=69;
step=69; % nonOverlap windows were used to crossvalidate the log Likelihood
TR=0.645; % this variable will not used for calculation, just for information
numVol=884;
numWin=floor((numVol-winSize)/step)+1;

% define file path
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
figDir=[analyDir,'fullBrainFC/figs/', covType, '/lambdaEstimate/'];
resultDir=[analyDir,'fullBrainFC/results/', covType, '/', session];

lambdaAtMaxLog=zeros(numSub,numWin,numSession);

for i=1:length(sessionList)
    session=char(sessionList{i});
    dataDir=[analyDir, '/data/645/all_10min/', covType, '/',session,'/'];
    
    parfor j=1:length(subList)
        subID=subList{j};
        
        disp(['Working on ',session, '_',subID, '.'])
        [logLikelihoodList, lambdaList]=lambdaEstimate(subID, winSize, step, numVol, dataDir);
        disp(['LogLikelihoodList and lambdaList for ', session, '_',subID, 'computed.'])
        
        for k=1:numWin
            logLikelihood=logLikelihoodList(k,:);
            lambda=lambdaList(k,:);
            Max=max(logLikelihood);
            lambdaAtMaxLog(j,k,i)=lambda(find(logLikelihood==Max));
        end
    end
end
meanLambdaAtMaxLogAcrossWin=mean(lambdaAtMaxLog,2)
save([resultDir,'/meanLambdaAtMaxLog_',session,'_',num2str(numSub),'subjects.mat'], 'meanLambdaAtMaxLogAcrossWin')
disp('lambdaAtMaxLog for all subjects and all sessions found and saved in one matrix.')


% plot the correlations between windows and subjects
% lambdaAtMaxLog2D=reshape(lambdaAtMaxLog,numSub,[]);
%
% corWinAllSession=corrcoef(lambdaAtMaxLog2D);
% figure(1)
% imagesc(corWinAllSession)
% title('corWinAllSession')
% colorbar('EastOutside')
% ylabel('LambdaAtMaxLog for Each Window')
% xlabel('LambdaAtMaxLog for Each Window')
% saveas(figure(1),[figDir, 'corWinAllSession.png']);
%
% corSub=corrcoef(lambdaAtMaxLog2D');
% figure(2)
% imagesc(corSub)
% title('corSub')
% colorbar('EastOutside')
% ylabel('LambdaAtMaxLog for Each Subject')
% xlabel('LambdaAtMaxLog for Each Subject')
% saveas(figure(2),[figDir, 'corSub.png']);
%
%
% % plot the distribution of lambda
% lambdaAtMaxLog1D=reshape(lambdaAtMaxLog,1,[]);
% figure(3)
% subplot(1,2,1)
% plot(lambdaAtMaxLog1D)
% ylabel('Lambda Value')
% xlabel('Number of Lambda')
% subplot(1,2,2)
% hist(lambdaAtMaxLog1D,0.05:0.01:0.2)
% ylabel('Frequency')
% xlabel('Lambda Value')
% saveas(figure(3),[figDir, 'Lambda Distribution.png']);

