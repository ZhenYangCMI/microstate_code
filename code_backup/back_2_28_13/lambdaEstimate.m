function [logLikelihoodList, lambdaList]=lambdaEstimate(session, subID)
% this function output the list for loglikelihood and lamda for a given sub
% and plot loglikelihood as a function of lamda
% input
% 1. session: string e.g. 'session1'
% 2. subID in cell format, e.g. subID={'0021002'};
% session='session1'
% subID={'0021002'}
winWidth=69;
step=69;
TR=0.645;
numVol=450;
numWin=floor((numVol-winWidth)/step)+1;


analyDir=(['/home/data/Projects/microstate/DPARSF_preprocessed/data/645/',session,'/']);
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
    disp (['Work on window ', num2str(q),' ......'])
    %         numVars = size(win(:,:,q), 2);
    %         sigma = cov(win(:,:,q));
    %
    %         tmp = max(max(abs(sigma)));
    %         lambdaListIn = (tmp ./ 10):(tmp ./ 10):tmp;
    clear tmp normWin
    tmp1=squeeze(win(:,:,q));
    normWin= (tmp1-repmat(mean(tmp1),size(tmp1,1),1))./repmat(std(tmp1),size(tmp1,1),1);
    normWin(isnan(normWin))=0;
    
    clear WListPerWin thetaListPerWin lambdaListPerWin errors
    %[WListPerWin, thetaListPerWin, lambdaListPerWin, errors] = GraphicalLassoPath(normWin);
    [WListPerWin, thetaListPerWin, lambdaListPerWin, errors] = ...
        GraphicalLassoPath(normWin,[0.01:0.01:1.0]);
    disp (['Graphical Lasso for window ', num2str(q),' done!'])
    
    
    %     % average covOtherWin across win
    %     if(k==1)
    %         meanCovOtherWin=covOtherWin;
    %     else
    %         temp=reshape(covOtherWin, [], numWin-1)';
    %         tempMean=mean(temp);
    %         meanCovOtherWin=reshape(tempMean,numROI,numROI);
    %     end
    %     disp(['Mean covariance matrix across other windows for window ', num2str(q), ' created!'])
    
    for m=1:length(lambdaListPerWin)
        
        % compute the average logLikelyhood for the rest win
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

% plot loglikelihood as a function of lambda and save the figs
% for q=1:numWin
%     figure(1)
%     subplot(2,3,q)
%     plot(lambdaList(q,:), logLikelihoodList(q,:))
% end
% 
% saveas(figure(1),[figDir,'lamdaEstimate_',session,'_sub_', char(subID), '.png'])
