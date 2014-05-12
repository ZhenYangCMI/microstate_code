function [winFullCor, winFullCorLasso, winPartialCorLasso]=FCWinCreationSingleSub(TR,session,sub)
% Output: this function output windows of full correlation or regularized covariance
% and inverse covariance matix for a single subject at multiple lambda
%Input:
%1.TRList: cell format, e.g. {'645','2500'};
%2.sessionlist: string, e.g. {'session1'}
%3.SubList: in cell format, e.g. subList={'0021002', '0021006'} or in a
%text file.
%4.corType: numerical, 1=full correlation; 2=regularized covariance (This
%is removed from the original function)


winWidth=69; % in TR: 69TR ~= 44s selected according to Allen's paper. 
step=3; % in TR
numVol=450;
numWin=floor((numVol-winWidth)/step)+1;
lambdaList=(0.08:0.01:0.14);
numLambda=length(lambdaList)

analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
resultDir=[analyDir,'results/',char(TR), '/',session,'/lambdaSensitivity/seedFCWin/', sub,'/'];

% create the convolved window
[finalWin]=winCreation(winWidth,0);
finalWin=finalWin';

disp(['Create full and partial windows for ',char(TR),' ', session,'.'])

disp (['Working on sub ', sub,' ......'])
subDir=[analyDir,'data/',char(TR),'/',session,'/',sub];
seedROISignals = load([subDir,'/ROISignals_seed_ROISignal.mat']);
TC1=seedROISignals.ROISignals;
numSeed=size(TC1,2);

ROIROISignals=load([subDir,'/ROISignals_atlas_ROISignal.mat']);
TC2=ROIROISignals.ROISignals;
numROI=size(TC2,2);

% concatenate the time series of seeds and ROIs
TC=[TC1,TC2];

% apply the sliding window to the time series
asize = size(TC);

numSub=1;
k=numSub;
% win(winWidth,asize(2),numWinAllSub);
for q=1+numWin*(k-1):numWin*k;
    for n=((q-1)-(k-1)*numWin)*step+1:winWidth+((q-1)-(k-1)*numWin)*step
        for m=1:asize(2)
            win(n-((q-1)-(k-1)*numWin)*step,m,q)=TC(n,m)*finalWin((n-((q-1)-(k-1)*numWin)*step),1);
        end
    end
end


disp (['Window applying done for all subjects. All windows of time series are saved in one matrix!'])
disp (['Compute the full and parital correlation for each window'])

% generate the full correlation for each window
for q=1:numWin
    corrWin(:,:,q)=corrcoef(win(:,:,q));
    for n=1:numSeed
        for m=1:numSeed
            dynFCBtwSeed(n,m,q)=corrWin(n,m,q);
            zDynFCBtwSeed(n,m,q)=0.5*log((1+dynFCBtwSeed(n,m,q))./(1-dynFCBtwSeed(n,m,q)));
        end
        for m=numSeed+1:numROI+numSeed
            dynFC(n,m-numSeed,q)=corrWin(n,m,q);
            zDynFC(n,m-numSeed,q)=0.5*log((1+dynFC(n,m-numSeed,q))./(1-dynFC(n,m-numSeed,q)));
        end
    end
end
disp('Full correlation between seeds and between seeds and ROIs are extracted for each window.')

winSeed1=reshape(zDynFC(1,:,:),numROI,[])';
winSeed2=reshape(zDynFC(2,:,:),numROI,[])';
winSeed3=reshape(zDynFC(3,:,:),numROI,[])';
winSeed4=reshape(zDynFC(4,:,:),numROI,[])';
winFullCor=vertcat(winSeed1, winSeed2, winSeed3, winSeed4);
disp ('Full correlation windows were concatenated as numWin per sub * numSub * numSeed.')
% save([resultDir,'winFullCor_sub', sub,'.mat'],'winFullCor')

% generate partial correlation for each window
for p=1:numLambda
    lambda=lambdaList(p)
    for q=1:numWin
        clear tmp normWin
        tmp1=squeeze(win(:,:,q));
        normWin= (tmp1-repmat(mean(tmp1),size(tmp1,1),1))./repmat(std(tmp1),size(tmp1,1),1);
        normWin(isnan(normWin))=0;
        [W, theta, lambdaOutput, errors]=GraphicalLassoPath(normWin,lambda);
        fullCorWin(:,:,q)=W;
        partialCorWin(:,:,q)=theta;
        disp(['Graphical Lasso for window',num2str(q),' done!'])
        for n=1:numSeed
            for m=1:numSeed
                partialDynFCBtwSeed(n,m,q)=partialCorWin(n,m,q);
                partialZDynFCBtwSeed(n,m,q)=0.5*log((1+partialDynFCBtwSeed(n,m,q))./(1-partialDynFCBtwSeed(n,m,q)));
            end
            for m=numSeed+1:numROI+numSeed
                partialDynFC(n,m-numSeed,q)=partialCorWin(n,m,q);
                partialZDynFC(n,m-numSeed,q)=0.5*log((1+partialDynFC(n,m-numSeed,q))./(1-partialDynFC(n,m-numSeed,q)));
                fullDynFC(n,m-numSeed,q)=fullCorWin(n,m,q);
                fullZDynFC(n,m-numSeed,q)=0.5*log((1+fullDynFC(n,m-numSeed,q))./(1-fullDynFC(n,m-numSeed,q)));
            end
        end
    end
    
    disp('Partial correlation between seeds and between seeds and ROIs are extracted for each window.')
    
    partWinSeed1=reshape(partialZDynFC(1,:,:),numROI,[])';
    partWinSeed2=reshape(partialZDynFC(2,:,:),numROI,[])';
    partWinSeed3=reshape(partialZDynFC(3,:,:),numROI,[])';
    partWinSeed4=reshape(partialZDynFC(4,:,:),numROI,[])';
    winPartialCorLasso=vertcat(partWinSeed1, partWinSeed2, partWinSeed3, partWinSeed4);
    
    fullWinSeed1=reshape(fullZDynFC(1,:,:),numROI,[])';
    fullWinSeed2=reshape(fullZDynFC(2,:,:),numROI,[])';
    fullWinSeed3=reshape(fullZDynFC(3,:,:),numROI,[])';
    fullWinSeed4=reshape(fullZDynFC(4,:,:),numROI,[])';
    winFullCorLasso=vertcat(fullWinSeed1, fullWinSeed2, fullWinSeed3, fullWinSeed4);
    
%     save([resultDir,'fullFCWin_sub', sub,'_lambda', num2str(lambda),'.mat'],'winFullCorLasso')
%     save([resultDir,'partrialFCWin_sub', sub,'_lambda', num2str(lambda),'.mat'],'winPartialCorLasso')
    
    disp ('Patial correlation windows were concatenated as numWin per sub * numSub * numSeed.')
end
end

