%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script estimate the optimum lambda for each subject based on the maximum log likelihood
%If two sessions were run, the lambdas for two sessions will be saved in one .mat file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

% Load the subjectList
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;
% subList={'0021002','0021018'};


% define the session
%sessionList={'session1','session2'};
sessionList={'session2'};

% define the window parameters
winSize=69; % in TRs
step=69; % in TRs, nonOverlap windows were used to crossvalidate the log Likelihood
TR=0.645; % in seconds. this variable will not used for calculation, just for information
numVol=884; % in TRs
numWin=floor((numVol-winSize)/step)+1;
ROISignalName=['ROISignals_TR2500k', numROI, '.mat'];

% define file path

resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/', covType, 'results/', session];

lambdaAtMaxLog=zeros(numSub,numWin,numSession);

for i=1:length(sessionList)
    session=char(sessionList{i});
    dataDir=['/home/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/', covType, '/',session,'/'];
    
    parfor j=1:length(subList)
        subID=subList{j};
        
        disp(['Working on ',session, '_',subID, '.'])
        [logLikelihoodList, lambdaList]=lambdaEstimate(subID, winSize, step, numVol, dataDir, ROISignalName);
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



