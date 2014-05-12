
%This script compute thess average sateionary FC across all
% subjects
clear
clc
close all

%TRList={'645','2500'};
sessionList={'session1','session2'};
subList=load('/home/data/Projects/microstate/NKITRT_SubID.mat');
subList=subList.SubID;
numSeed=4;
numROI=156;

TRList={'645'};
%sessionList={'session1'};
% subList={'0021002'};

numSub=length(subList);
numSession=length(sessionList);
numTR=length(TRList);
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];

for i=1:numTR
    TR=char(TRList{i});
    for j=1:numSession
        session=char(sessionList{j});
        resultDir=[analyDir,'results/',TR, '/all_10min/',session,'/'];
        correl=zeros(numSeed+numROI,numSeed+numROI,numSub);
        stat_FC_seed=zeros(numSeed, numSeed, numSub);
        stat_FC=zeros(numSeed, numROI, numSub);
        for k=1:numSub
            sub=subList{k};
            disp (['Working on sub ', char(sub),' ......'])
            subDir=[analyDir,'data/',TR,'/all_10min/',session,'/',char(sub)];
            seedROISignals = load([subDir,'/ROISignals_seedROISignal.mat']);
            TC1=seedROISignals.ROISignals;
            %numSeed=size(TC1,2);
            
            ROIROISignals=load([subDir,'/ROISignals_atlasROISignal.mat']);
            TC2=ROIROISignals.ROISignals;
            %numROI=size(TC2,2);
            
            % concatenate the time series of seeds and ROIs
            TC=[TC1,TC2];
            
            % apply the sliding window to the time series
            asize = size(TC);
            
            %Compute the full correlation matrix
            correl(:,:,k)=corrcoef(TC);
            
            % Extract the stationary FC between all pairs of seeds
            for n=1:numSeed
                for m=1:numSeed
                    stat_FC_seed(n,m,k)=correl(n,m,k);
                    z_stationary_FC_seed(n,m,k)=0.5*log((1+stat_FC_seed(n,m,k))./(1-stat_FC_seed(n,m,k)));
                end
            end
            
            % Extract the stationary FC between seeds and ROIs
            for n=1:numSeed
                for m=numSeed+1:numROI+numSeed
                    stat_FC(n,m-numSeed,k)=correl(n,m,k);
                    z_stationary_FC(n,m-numSeed,k)=0.5*log((1+stat_FC(n,m-numSeed,k))./(1-stat_FC(n,m-numSeed,k)));
                end
            end
            
        end
        
        % average across all subjects
        tmp1=reshape(z_stationary_FC_seed,[],k)';
        avg_tmp1=mean(tmp1);
        avg_z_stationary_FC_seed=reshape(avg_tmp1,numSeed,numSeed);
        save([resultDir,'/avg_z_stationary_FC_seed.mat'],'avg_z_stationary_FC_seed');
        save([resultDir,'/z_stationary_FC_seed.mat'],'z_stationary_FC_seed');
        
        tmp2=reshape(z_stationary_FC,[],k)';
        avg_tmp2=mean(tmp2);
        avg_z_stationary_FC=reshape(avg_tmp2,numSeed,numROI);
        save([resultDir,'/avg_z_stationary_FC.mat'],'avg_z_stationary_FC');
        save([resultDir,'/z_stationary_FC.mat'],'z_stationary_FC');
    end
end