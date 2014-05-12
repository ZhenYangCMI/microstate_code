function [winFullCorSeed, winFullCorLassoSeed, winPartialCorLassoSeed]=FCWinCreationBtwSeed(TR,session,subList,lambdaList, covType)
% Output: this function output windows of full correlation or regularized covariance
% with all subjects and all seeds concatenated.
%Input:
%1.TR: string  e.g. '645','2500';
%2.session: string e.g. 'session1'
%3.SubList: in cell format, e.g. subList={'0021002', '0021006'} or in a
%text file.
%4.lambdaList: numerical e.g. [0.11 0.12]

% session='session1'
% TR='645';
% subList={'0021002','0021006'};
% lambdaList=[0.11, 0.12];
% covType='GSR';

% in TR: 69TR ~= 44s selected according to Allen's paper and was used in the main analysis. 22s ~=34TR 88~=136TR were also tested 
% winWidth=69; 

winWidth=69;
step=3; % in TR
numVol=884;
numWin=floor((numVol-winWidth)/step)+1;
numSub=length(subList);
numWinAllSub=numWin*numSub;
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];

% create the convolved window
[finalWin]=winCreation(winWidth,0);
finalWin=finalWin';

disp(['Create full and partial correlation windows for ',TR,' ', session,'.'])

for k=1:numSub
    sub=subList{k};
    disp (['Working on sub ', char(sub),' ......'])
    subDir=[analyDir,'data/',TR,'/all_10min/',covType, '/',session,'/',char(sub)];
    seedROISignals = load([subDir,'/ROISignals_seedROISignal.mat']);
    TC1=seedROISignals.ROISignals;
    numSeed=size(TC1,2);
    
    ROIROISignals=load([subDir,'/ROISignals_atlasROISignal.mat']);
    TC2=ROIROISignals.ROISignals;
    numROI=size(TC2,2);
    
    % concatenate the time series of seeds and ROIs
    TC=[TC1,TC2];
    
    % apply the sliding window to the time series
    asize = size(TC);
    
    % win(winWidth,asize(2),numWinAllSub);
    for q=1+numWin*(k-1):numWin*k;
        for n=((q-1)-(k-1)*numWin)*step+1:winWidth+((q-1)-(k-1)*numWin)*step
            for m=1:asize(2)
                win(n-((q-1)-(k-1)*numWin)*step,m,q)=TC(n,m)*finalWin((n-((q-1)-(k-1)*numWin)*step),1);
            end
        end
    end
end

disp (['Window applying done for all subjects. All windows of time series are saved in one matrix!'])
disp (['Compute the full and parital correlation for each window'])

% generate the full correlation for each window
for q=1:numWinAllSub
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

winAllSubSeed1=reshape(dynFCBtwSeed(1,:,:),numSeed,[])';
winAllSubSeed2=reshape(dynFCBtwSeed(2,:,:),numSeed,[])';
winAllSubSeed3=reshape(dynFCBtwSeed(3,:,:),numSeed,[])';
winAllSubSeed4=reshape(dynFCBtwSeed(4,:,:),numSeed,[])';
winFullCorSeed=vertcat(winAllSubSeed1, winAllSubSeed2, winAllSubSeed3, winAllSubSeed4);

disp ('Full correlation windows were concatenated as numWin/sub * numSub * numSeed.')

% generate partial correlation for each window
for q=1:numWinAllSub
    clear tmp normWin
    tmp1=squeeze(win(:,:,q));
    normWin= (tmp1-repmat(mean(tmp1),size(tmp1,1),1))./repmat(std(tmp1),size(tmp1,1),1);
    normWin(isnan(normWin))=0;
    lambdaIndx=ceil(q/numWin);
    lambda=lambdaList(lambdaIndx)
    [W, theta, lambdaOut, errors]=GraphicalLassoPath(normWin,lambda);
    fullCorWin(:,:,q)=W;
    partialCorWin(:,:,q)=theta;
    disp(['Graphical Lasso for window',num2str(q),' done!'])
    for n=1:numSeed
        for m=1:numSeed
            partialDynFCBtwSeed(n,m,q)=partialCorWin(n,m,q);
            partialZDynFCBtwSeed(n,m,q)=0.5*log((1+partialDynFCBtwSeed(n,m,q))./(1-partialDynFCBtwSeed(n,m,q)));
            
            fullDynFCBtwSeed(n,m,q)=fullCorWin(n,m,q);
            fullZDynFCBtwSeed(n,m,q)=0.5*log((1+fullDynFCBtwSeed(n,m,q))./(1-fullDynFCBtwSeed(n,m,q)));
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

partWinAllSubSeed1=reshape(partialDynFCBtwSeed(1,:,:),numSeed,[])';
partWinAllSubSeed2=reshape(partialDynFCBtwSeed(2,:,:),numSeed,[])';
partWinAllSubSeed3=reshape(partialDynFCBtwSeed(3,:,:),numSeed,[])';
partWinAllSubSeed4=reshape(partialDynFCBtwSeed(4,:,:),numSeed,[])';
winPartialCorLassoSeed=vertcat(partWinAllSubSeed1, partWinAllSubSeed2, partWinAllSubSeed3, partWinAllSubSeed4);

fullWinSeed1=reshape(fullDynFCBtwSeed(1,:,:),numSeed,[])';
fullWinSeed2=reshape(fullDynFCBtwSeed(2,:,:),numSeed,[])';
fullWinSeed3=reshape(fullDynFCBtwSeed(3,:,:),numSeed,[])';
fullWinSeed4=reshape(fullDynFCBtwSeed(4,:,:),numSeed,[])';
winFullCorLassoSeed=vertcat(fullWinSeed1, fullWinSeed2, fullWinSeed3, fullWinSeed4);

disp ('Patial correlation windows were concatenated as numWin/sub * numSub * numSeed.')
end




