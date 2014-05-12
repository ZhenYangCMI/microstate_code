function [zFullCorrWin, zfullCorrWinLasso, zpatialCorrWinLasso]=FCWinCreation(dataDir,subList,lambdaList, winSize, step, numVol, ROISignalName)
% Output: this function output windows of full correlation or regularized covariance
% with all subjects and all seeds concatenated.
%Input:
%1.dataDir: string path to the sub TS
%2.SubList: in cell format, e.g. subList={'0021002', '0021006'} or in a
%text file.
%3.lambdaList: numerical e.g. [0.11 0.12]
%4. winSize: width of the wind size
%5. ROISignalName: the name of the ROI signal, .mat file

% e.g.
% dataDir=dataDir=[analyDir, '/data/645/all_10min/', covType, '/',session,'/'];
% subList={'0021002'};
% lambdaList=[0.11];
% winSize=69;


% in TR: 69TR ~= 44s selected according to Allen's paper and was used in the main analysis. 22s ~=34TR 88~=136TR were also tested
% winSize=69;

numWin=floor((numVol-winSize)/step)+1;
numSub=length(subList);
numWinAllSub=numWin*numSub;

% create the convolved window
[finalWin]=winCreation(winSize,0);
finalWin=finalWin';

% segment the TS into TS windows

for k=1:numSub
    sub=subList{k};
    disp (['Working on sub ', char(sub),' ......'])
    subDir=[dataDir,char(sub)];
    ROISignals = load([subDir,'/', ROISignalName]);
    TC=ROISignals.ROISignals;
    
    % apply the sliding window to the time series
    asize = size(TC);
    
    % win(winSize,asize(2),numWinAllSub);
    for q=1+numWin*(k-1):numWin*k;
        for n=((q-1)-(k-1)*numWin)*step+1:winSize+((q-1)-(k-1)*numWin)*step
            for m=1:asize(2)
                win(n-((q-1)-(k-1)*numWin)*step,m,q)=TC(n,m)*finalWin((n-((q-1)-(k-1)*numWin)*step),1);
            end
        end
    end
end
disp (['Window applying done for all subjects. All windows of time series are saved in one matrix!'])


% generate full and partial correlation for each window
disp (['Compute the full and parital correlation for each window'])
for q=1:numWinAllSub
    % compute the full correlation: Pearson's r
    fullCorrWin(:,:,q)=corrcoef(win(:,:,q));
    % compute the full and partial correlation: Lasso
    clear tmp normWin
    tmp1=squeeze(win(:,:,q));
    normWin= (tmp1-repmat(mean(tmp1),size(tmp1,1),1))./repmat(std(tmp1),size(tmp1,1),1);
    normWin(isnan(normWin))=0;
    lambdaIndx=ceil(q/numWin);
    lambda=lambdaList(lambdaIndx)
    [W, theta, lambdaOut, errors]=GraphicalLassoPath(normWin,lambda);
    fullCorrWinLasso(:,:,q)=W;
    partialCorrWinLasso(:,:,q)=theta;
    disp(['Graphical Lasso for window',num2str(q),' done!'])
end

% Fisher z tranform the correlations
zFullCorrWin=0.5*log((1+fullCorrWin)./(1-fullCorrWin));
zfullCorrWinLasso=0.5*log((1+fullCorrWinLasso)./(1-fullCorrWinLasso));
zpatialCorrWinLasso=0.5*log((1+partialCorrWinLasso)./(1-partialCorrWinLasso));
disp('Full and partial correlation are computed for each window.')

end




