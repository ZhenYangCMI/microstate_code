function [winAllSubAllSeed]=FCWinCreation(TRList,sessionList,subList, corType)
% Output: this function output windows of full correlation or regularized covariance
% with all subjects and all seeds concatenated.
%Input:
%1.TRList: cell format, e.g. {'645','2500'};
%2.sessionlist: string, e.g. {'session1'}
%3.SubList: in cell format, e.g. subList={'0021002', '0021006'} or in a
%text file.
%4.corType: numerical, 1=full correlation; 2=regularized covariance


winWidth=69;
step=3;
numVol=450;
numWin=floor((numVol-winWidth)/step)+1;
numSub=length(subList);
numSession=length(sessionList);
numTR=length(TRList);
numWinAllSub=numWin*numSub;

%analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
analyDir=['/Users/zhenyang/Documents/microstate/'];

dataDir=([analyDir,'data']);
figDir=([analyDir,'fig']);
resultDir=([analyDir,'results']);
maskDir=([analyDir,'mask/']);


for i=1:numTR
    TR=TRList{i};
    TRDir=[dataDir,'/',TR];
    if ~exist([resultDir,'/',TR], 'dir')
        mkdir(resultDir,TR)
    end
    resultTRDir=[resultDir,'/',TR];
    
    for j=1:numSession
        session=sessionList{j};
        sessionDir=[TRDir,'/',session];
        if ~exist([resultTRDir,'/',session], 'dir')
            mkdir(resultTRDir, session)
        end
        resultSesDir=[resultTRDir,'/',session];
        
        for k=1:numSub
            subDir=[sessionDir,'/', char(subList{k})];
            disp (['Working on sub ', char(subList{k}),' ......'])
            
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
            [finalWin]=winCreation(winWidth,0);
            finalWin=finalWin';
            
            % win(W_width,asize(2),N_win_tot);
            for q=1+numWin*(k-1):numWin*k;
                for n=((q-1)-(k-1)*numWin)*step+1:winWidth+((q-1)-(k-1)*numWin)*step
                    for m=1:asize(2)
                        win(n-((q-1)-(k-1)*numWin)*step,m,q)=TC(n,m)*finalWin((n-((q-1)-(k-1)*numWin)*step),1);
                    end
                end
            end
        end
        
        
        if (corType==1)
            % generate the full correlation matrix for each win and extract the dyn FC
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
            save([resultSesDir,'/zDynFCFullCorr.mat'],'zDynFC');
            disp('Full correlation between seeds and between seeds and ROIs are extracted for each window.')
            
            winAllSubSeed1=reshape(zDynFC(1,:,:),numROI,[])';
            winAllSubSeed2=reshape(zDynFC(2,:,:),numROI,[])';
            winAllSubSeed3=reshape(zDynFC(3,:,:),numROI,[])';
            winAllSubSeed4=reshape(zDynFC(4,:,:),numROI,[])';
            winAllSubAllSeed=vertcat(winAllSubSeed1, winAllSubSeed2, winAllSubSeed3, winAllSubSeed4);
            disp ('Windows were concatenated as numWin/sub * numSub * numSeed.')
            save([resultSesDir,'/winAllSubAllSeed_fullCor_',char(TR),'_',char(session),'.mat'],'winAllSubAllSeed');
        
        else
            for q=1:numWinAllSub
                partialCorWin(:,:,q)=GraphicLasso(squeeze(win(:,:,q),0.1);
                
                for n=1:numSeed
                    for m=1:numSeed
                        partialDynFCBtwSeed(n,m,q)=partialCorWin(n,m,q);
                        partialZDynFCBtwSeed(n,m,q)=0.5*log((1+partialDynFCBtwSeed(n,m,q))./(1-partialDynFCBtwSeed(n,m,q)));
                    end
                    for m=numSeed+1:numROI+numSeed
                        partialDynFC(n,m-numSeed,q)=partialCorWin(n,m,q);
                        partialZDynFC(n,m-numSeed,q)=0.5*log((1+partialDynFC(n,m-numSeed,q))./(1-paritalDynFC(n,m-numSeed,q)));
                    end
                end
            end
            save([resultSesDir,'/zDynFCPartialCorr.mat'],'partialZDynFC');
            disp('Partial correlation between seeds and between seeds and ROIs are extracted for each window.')
            
            winAllSubSeed1=reshape(partialZDynFC(1,:,:),numROI,[])';
            winAllSubSeed2=reshape(partialZDynFC(2,:,:),numROI,[])';
            winAllSubSeed3=reshape(partialZDynFC(3,:,:),numROI,[])';
            winAllSubSeed4=reshape(partialZDynFC(4,:,:),numROI,[])';
            winAllSubAllSeed=vertcat(winAllSubSeed1, winAllSubSeed2, winAllSubSeed3, winAllSubSeed4);
            disp ('Windows were concatenated as numWin/sub * numSub * numSeed.')
            save([resultSesDir,'/winAllSubAllSeed_partialCor_',char(TR),'_',char(session),'.mat'],'winAllSubAllSeed');
        end
    end
end

