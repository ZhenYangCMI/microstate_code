function [logLikelihoodList, lambdaList]=lambdaEstimate(subID, winSize, step, numVol, dataDir, ROISignalName)
% this function output the loglikelihood list and the lamda list for a given sub and all possible windows
% and plot loglikelihood as a function of lamda
% input
% 1. subID in cell format, e.g. subID={'0021002'};
% 2. winSize, numeric, in TRs
% 3. step, numeric, in TRs
% 4. numVol: number of volumes
% 5. dataDir: string, the data Dir, one level up the sub dir

%% below are example for testing the function
% clear
% clc
% 
% subID={'0021002'}
% winSize=69;
% step=69;
% TR=0.645;
% numVol=884;
% dataDir=(['/home/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/GSR/session1/']);
%%
%%

numWin=floor((numVol-winSize)/step)+1;

subDir=[dataDir,char(subID)];

% read in the TS for this sub
ROISignals = load([subDir,'/', ROISignalName]);
TC=ROISignals.ROISignals;
numROI=size(TC,2);

% create the sliding window (gaussian convolve rectangular)
[finalWin]=winCreation(winSize,0);
finalWin=finalWin';

% apply the sliding window to the time series
win=zeros(winSize, numROI, numWin);
for q=1:numWin
    for n=(q-1)*step+1:winSize+(q-1)*step
        for m=1:size(TC,2)
            win(n-(q-1)*step,m,q)=TC(n,m)*finalWin(n-(q-1)*step,1);
        end
    end
end

% for each win, generate the thetaList and lamdaList
for q=1:numWin
    disp (['Work on window ', num2str(q),' ......'])
    
    clear tmp normWin
    tmp1=squeeze(win(:,:,q));
    % Z standardize the TS of each window before run gLasso. It will not work if don't standardize.
    normWin= (tmp1-repmat(mean(tmp1),size(tmp1,1),1))./repmat(std(tmp1),size(tmp1,1),1);
    normWin(isnan(normWin))=0;
    
    clear WListPerWin thetaListPerWin lambdaListPerWin errors
    % here PerWin means the num of nonoverlapping windows used to do crossvalidation.
    [WListPerWin, thetaListPerWin, lambdaListPerWin, errors] = ...
        GraphicalLassoPath(normWin,[0.01:0.01:1.0]);
    disp (['Graphical Lasso for window ', num2str(q),' done!'])
    
    % this is average the cov across windows then compute the log likelihood, the other way seems more appropriate
    %     % average covOtherWin across win
    %     if(k==1)
    %         meanCovOtherWin=covOtherWin;
    %     else
    %         temp=reshape(covOtherWin, [], numWin-1)';
    %         tempMean=mean(temp);
    %         meanCovOtherWin=reshape(tempMean,numROI,numROI);
    %     end
    %     disp(['Mean covariance matrix across other windows for window ', num2str(q), ' created!'])
    
    % This is compute the log likelihood first, then average the log likelihood
    for m=1:length(lambdaListPerWin)
        
        % compute the cov matrix for the validation window
        k=0;
        for j=1:numWin
            if (j~=q)
                k=k+1;
                clear tmp2 normOtherWin
                tmp2=squeeze(win(:,:,j));
                normOtherWin= (tmp2-repmat(mean(tmp2),size(tmp2,1),1))./repmat(std(tmp2),size(tmp2,1),1);
                normOtherWin(isnan(normOtherWin))=0;
                covOtherWin(:,:,k)=cov(normOtherWin);
            end
        end
        
        % compute the log likelihood, k equals to the number of windows other than the target window used to compute the theta
        clear perwin
        thetaPerWin=thetaListPerWin(:,:,m);
        
        for n=1:k
            covPerWin=covOtherWin(:,:,n);
            logLikelihood(n,m)=log(det(thetaPerWin))...
                -trace(covPerWin*thetaPerWin)...
                -lambdaListPerWin(m)*norm(thetaPerWin,1);
        end
        meanLogLikelihood=mean(logLikelihood);
        logLikelihoodList(q,m)=meanLogLikelihood(1,m);
        lambdaList(q,m)=lambdaListPerWin(m);
    end
    disp(['loglikelihood for win ', num2str(q), ' computed!'])
end
disp ('loglikelihood and lamda list created!')
end


