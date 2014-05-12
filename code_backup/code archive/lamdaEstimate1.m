function [logLikelihoodList, lamdaList]=lamdaEstimate(session, subID)
% this function output the list for loglikelihood and lamda for a given sub
% and plot loglikelihood as a function of lamda
% input
% 1. session: string e.g. 'session1'
% 2. subID in cell format, e.g. subID={'0021002'};
session='session1'
subID={'7055197'}
winWidth=69;
step=69;
TR=0.645;
numVol=450;
numWin=floor((numVol-winWidth)/step)+1;


analyDir=(['/home/data/Projects/microstate/DPARSF_preprocessed/data/645/first_5min/',session,'/']);
subDir=[analyDir,char(subID)];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/'];


seedROISignals = load([subDir,'/ROISignals_seed_ROISignal.mat']);
TC1=seedROISignals.ROISignals;

ROIROISignals=load([subDir,'/ROISignals_atlas_ROISignal.mat']);
TC2=ROIROISignals.ROISignals;

% concatenate the time series of seeds and ROIs
TC=[TC1,TC2];

% apply the sliding window to the time series
numROI=size(TC,2);

% apply the sliding window to the time series
[finalWin]=winCreation(winWidth,0);
finalWin=finalWin';


for q=1:numWin
    for n=(q-1)*step+1:winWidth+(q-1)*step
        for m=1:size(TC,2)
            win(n-(q-1)*step,m,q)=TC(n,m)*finalWin(n-(q-1)*step,1);
        end
    end
end

% for each win, generate the thetaList and lamdaList
%for q=1:numWin
for q=1:numWin
%         numVars = size(win(:,:,q), 2);
%         sigma = cov(win(:,:,q));
%     
%         tmp = max(max(abs(sigma)));
%         lambdaListIn = (tmp ./ 10):(tmp ./ 10):tmp;
    
     tmp=win(:,:,q);
    normWin= (tmp-repmat(mean(tmp),size(tmp,1),1))./repmat(std(tmp),size(tmp,1),1);
    normWin(isnan(normWin))=0;

    disp (['Work on window ', num2str(q),' ......'])
    clear WListPerWin thetaListPerWin lambdaListPerWin errors
    
    psize=size(win);
    for i=1:psize(1,1)
        for j=1:psize(1,2)
            qwin(i,j)=win(i,j,q);
        end
    end
   
    
    [WListPerWin, thetaListPerWin, lambdaListPerWin, errors] = GraphicalLassoPath(qwin, lambdaListIn,0,0,1,1e-4,1e4);
    
    % compute the average covariance matrix for the rest win
    k=0;
    for j=1:numWin
        if (j~=q)
            k=k+1;
            covOtherWin(:,:,k)=cov(squeeze(win(:,:,j)));
        end
    end
    
    % average covOtherWin across win
    if(k==1)
        meanCovOtherWin=covOtherWin;
    else
        temp=reshape(covOtherWin, [], numWin-1)';
        tempMean=mean(temp);
        meanCovOtherWin=reshape(tempMean,numROI,numROI);
    end
    disp(['Mean covariance matrix across other windows for window ', num2str(q), ' created!'])
    
    for m=1:length(lambdaListPerWin)
        qsize=size(thetaListPerWin);
        for i=1:qsize(1,1)
            for j=1:qsize(1,2)
                perwin(i,j)=thetaListPerWin(i,j,m);
            end
        end
        logLikelihoodList(q,m)=log(det(perwin))...
            -trace(meanCovOtherWin*perwin)...
            -lambdaListPerWin(m)*norm(perwin,1);
        lambdaList(q,m)=lambdaListPerWin(m);
    end
    disp(['loglikelihood for win ', num2str(q), ' computed!'])
end
disp ('loglikelihood and lamda list created!')

for q=1:numWin
    figure(1)
    subplot(2,3,q)
    plot(lambdaList(q,:), logLikelihoodList(q,:))
end

saveas(figure(1),[figDir,'lamdaEstimate_session',session,'_sub_', char(subID), '.png'])
