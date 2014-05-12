function [logLikelihoodList, lamdaList]=lamdaEstimate(subID)
% this function output the list for loglikelihood and lamda for a given sub
% input subID in cell format, e.g. subID={'0021002'};

winWidth=69;
step=69;
TR=0.645;
numVol=450;
numWin=floor((numVol-winWidth)/step)+1;


analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/data/645/first_5min/session1/'];
subDir=[analyDir,char(subID)];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/'];
disp (['Working on sub ', char(subID), ' ......'])

seedROISignals = load([subDir,'/ROISignals_seed_ROISignal.mat']);
TC1=seedROISignals.ROISignals;

ROIROISignals=load([subDir,'/ROISignals_atlas_ROISignal.mat']);
TC2=ROIROISignals.ROISignals;

% concatenate the time series of seeds and ROIs
TC=[TC1,TC2];

% apply the sliding window to the time series
numROI=size(TC,2);

% apply the sliding window to the time series
[finalWin]=winCreation(winWidth);
finalWin=finalWin';

for q=1:numWin
    for n=(q-1)*step+1:winWidth+(q-1)*step
        for m=1:size(TC,2)
            win(n-(q-1)*step,m,q)=TC(n,m)*finalWin(n-(q-1)*step,1);
        end
    end
end

% for each win, generate the thetaList and lamdaList
for i=1:numWin
    
    [WListPerWin, thetaListPerWin, lamdaListPerWin, errors]=GraphicalLassoPath(squeeze(win(:,:,i)));
    
    % compute the average covariance matrix for the rest win
    k=0
    for j=1:numWin
        if (j~=i)
            k=k+1;
            covOtherWin(:,:,k)=cov(squeeze(win(:,:,j)))
            
            % average covOtherWin across win
            temp=reshape(covOtherWin, [], numWin-1)';
            tempMean=mean(temp);
            meanCovOtherWin=reshape(tempMean,numROI,numROI);
        end
    end
    for m=1:length(lamdaListPerWin)
        logLikelihood(i,m)=log(det(thetaListPerWin(:,:,m)))...
            -trace(meanCovOtherWin*thetaListPerWin(:,:,m))...
            -lamdaListPerWin(m)*(thetaListPerWin(:,:,m))
        lamdaList(i,m)=lamdaListPerWin(m)
    end
end
    
for i=1:numWin
    figure(1)
    subplot(2,3,i)
    plot(lamdaList(i,:), logLikelihood(i,:))
end

saveas(figure(1),[figDir,'lamdaEstimate_sub_', char(subID), '.png'])
        
        
        
        
        
        
        
        
        
        
        
        
        %