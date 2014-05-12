%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script estimate the optimum lambda for each subject based on the maximum log likelihood
%If two sessions were run, the lambdas for two sessions will be saved in one .mat file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

covType='GSR';

% Load the subjectList
%subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
%subList=subList.SubID;
subList={'3808535','8574662'};
numSub=length(subList);

% define the session
%sessionList={'session1','session2'};
session='session1';
numROI='0200'
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];

% define the window parameters
winSize=69;
step=3; % in TR
TR=0.645;
numVol=442; % 5min data
numWin=floor((numVol-winSize)/step)+1;
ROISignalName=['ROISignals_TR645k', numROI, '_5min.mat'];
dataDir=[analyDir, '/fullBrainFC/testing/TR645_5min/'];
figDir=[analyDir,'fullBrainFC/testing/figs/TR645_5min/'];
resultDir=[analyDir,'fullBrainFC/testing/results/TR645_5min/'];

% below is for TR 2500
% winSize=18;
% step=1; % nonOverlap windows were used to crossvalidate the log Likelihood
% TR=2.5; % this variable will not used for calculation, just for information
% numVol=116;
% numWin=floor((numVol-winSize)/step)+1;
% ROISignalName=['ROISignals_TR2500k', numROI, '.mat'];

lambdaAtMaxLog=zeros(numSub,numWin);

for j=1:length(subList)
    subID=subList{j};
    
    disp(['Working on ',session, '_',subID, '.'])
    [logLikelihoodList, lambdaList]=lambdaEstimate(subID, winSize, step, numVol, dataDir, ROISignalName);
    disp(['LogLikelihoodList and lambdaList for ', session, '_',subID, 'computed.'])
    
    for k=1:numWin
        logLikelihood=logLikelihoodList(k,:);
        lambda=lambdaList(k,:);
        Max=max(logLikelihood);
        lambdaAtMaxLog(j,k)=lambda(find(logLikelihood==Max));
    end
end

meanLambdaAtMaxLogAcrossWin=mean(lambdaAtMaxLog,2)
save([resultDir,'/meanLambdaAtMaxLog_',numROI, '_', session,'_',num2str(numSub),'subjects_TR', num2str(TR), '.mat'], 'meanLambdaAtMaxLogAcrossWin')
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

