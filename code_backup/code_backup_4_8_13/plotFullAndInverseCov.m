

clear
clc

session='session1';
subID='8574662';
lambdaList=[0.08:0.01:0.14];
lambda=0.12;
optimalLambdaIndx=find(lambdaList==lambda);

maskDir=['/home/data/Projects/microstate/DPARSF_preprocessed/mask/'];
resultsDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',session,'/lambdaSensitivity/seedFCWin/',subID,'/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',session,'/lambdaSensitivity/',subID,'/'];

tmp1=load([resultsDir, 'partrialFCWin_sub',subID,'_lambda',num2str(lambda),'.mat']);
tmp2=load([resultsDir, 'winFullCor_sub',subID,'.mat']);
tmp3=load([resultsDir, 'correlEachSeedType2.mat']);

 winPartialCorLasso=tmp1.winPartialCorLasso;
 winFullCor=tmp2.winFullCor;
 correlEachSeedType2=tmp3.correlEachSeedType2;
 
 % extract the corrcoef for all windows on optimaLambda (lambdaIndx=5) for seed 2
 
 correlOptimaLambda=squeeze(correlEachSeedType2(:,5,2));
 minIndx=find(correlOptimaLambda==max(correlOptimaLambda));
 maxIndx=find(correlOptimaLambda==min(correlOptimaLambda));
 
 winFullCor=winFullCor(129:256,:);
 meanWinFullCor=mean(winFullCor);
 maxWinFullCor=winFullCor(maxIndx,:);
 minWinFullCor=winFullCor(minIndx,:);
 
 inverseCov=winPartialCorLasso(129:256,:);
 winPartialCor=-inverseCov;
 meanWinPartialCor=mean(winPartialCor);
 maxWinPartialCor=winPartialCor(maxIndx,:);
 minWinPartialCor=winPartialCor(minIndx,:);
 
 disp('Create brain map of windows')
      
        [Outdata,VoxDim,Header]=rest_readfile([maskDir,'final_reduced.nii']);
        [nDim1 nDim2 nDim3]=size(Outdata);
        temp=unique(Outdata);
        ROIIndx=temp(find(temp~=0));
        numROI=length(ROIIndx);
        
        stateMap1=Outdata;
        stateMap2=Outdata;
        stateMap3=Outdata;
        stateMap4=Outdata;
        stateMap5=Outdata;
        stateMap6=Outdata;
        for m=1:numROI
            stateMap1(find(Outdata==ROIIndx(m)))=meanWinFullCor(m);
            stateMap2(find(Outdata==ROIIndx(m)))=maxWinFullCor(m);
            stateMap3(find(Outdata==ROIIndx(m)))=minWinFullCor(m);
            stateMap4(find(Outdata==ROIIndx(m)))=meanWinPartialCor(m);
            stateMap5(find(Outdata==ROIIndx(m)))=meanWinPartialCor(m);
            stateMap6(find(Outdata==ROIIndx(m)))=meanWinPartialCor(m);
        end
        
        
        Header.pinfo = [1;0;0];
        Header.dt    =[16,0];
        rest_WriteNiftiImage(stateMap1,Header,[figDir,'meanWinFullCor.nii']);
        rest_WriteNiftiImage(stateMap2,Header,[figDir,'maxWinFullCor.nii']);
        rest_WriteNiftiImage(stateMap3,Header,[figDir,'minWinFullCor.nii']);
        rest_WriteNiftiImage(stateMap4,Header,[figDir,'meanWinPartialCor.nii']);
        rest_WriteNiftiImage(stateMap5,Header,[figDir,'maxWinPartialCor.nii']);
        rest_WriteNiftiImage(stateMap6,Header,[figDir,'minWinPartialCor.nii']);
        
        
        
        
        
        
        
        
        
 