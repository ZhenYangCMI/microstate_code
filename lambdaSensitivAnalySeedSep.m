function [correlEachSeedType3, distEachSeedType3,mutualInfoEachSeedType3]=lambdaSensitivAnalySeedSep(TR,session, sub, optimalLambda)
%the following function estimate the sensitivity of lambda to the FC
%windows using 3 metric: correlation, euclidean distance, and mutual
%information. There are three types of comparisons across 7 lambdas. Type1: compare the similarity between the full correlation matrix generated using Lasso
%(i.e. covariance matrix computed on normalized TS of that window; Type2:
%compare the full correlation matrix generated using matlab corrcoef function and
%the inverse covariance matrix generated using GLasso on normalized TS of
%that window (i.e. partial correlation); Type 3; compare the parital
%correlation matrix at differen lambda with the partial correlation matrix
%computed at the lambdaAtMaxLog
         

%1 TR: cell, e.g. TR={'645'}; 2mutualInfoEachSeedType1 session: string, e.g 'session1'; 3.
% sub:string, e.g. '8574662'; 4. optimalLambda: numerical

% clear
% clc
% 
% sub='8574662';
% session='session1'
% TR={'645'};
% optimalLambda=0.12;


numWin=128;
lambdaList=(0.08:0.01:0.14);
numLambda=length(lambdaList)
numSeed=4;

analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
resultDir=[analyDir,'results/',char(TR), '/',session,'/lambdaSensitivity/seedFCWin/', sub,'/'];

% evaluate the sensitivity for each seed
tmp=load([resultDir,'partrialFCWin_sub', sub,'_lambda', num2str(optimalLambda),'.mat'])
optimalWinPartialCorLasso =tmp.winPartialCorLasso;

for i=1:numLambda
    lambda=lambdaList(i)
    tmp1=load([resultDir,'winFullCor_sub', sub,'.mat'])
    tmp2=load([resultDir,'fullFCWin_sub', sub,'_lambda', num2str(lambda),'.mat'])
    tmp3=load([resultDir,'partrialFCWin_sub', sub,'_lambda', num2str(lambda),'.mat'])
    
    winFullCor=tmp1.winFullCor;
    winFullCorLasso=tmp2.winFullCorLasso;
    winPartialCorLasso=tmp3.winPartialCorLasso;
    
    % estimate the similarity between the full correlation matrix and the covariance matrix, the full correlation matrix and the inversed covariance matrix
    
    for k=1:numSeed
        disp(['Working on seed ', num2str(k)])
        winFullCorSeed=winFullCor(1+(k-1)*numWin:k*numWin,:);
        winFullCorLassoSeed= winFullCorLasso(1+(k-1)*numWin:k*numWin, :);
        winPartialCorLassoSeed=winPartialCorLasso(1+(k-1)*numWin:k*numWin,:);
        optimalWinPartialCorLassoSeed=optimalWinPartialCorLasso(1+(k-1)*numWin:k*numWin,:);
        for j=1:numWin
            % use correlation to estimate similarity
            disp(['Computing correlations for window ',num2str(j)])
            correl1=corrcoef(winFullCorSeed(j,:),winFullCorLassoSeed(j,:));
            correl2=corrcoef(winFullCorSeed(j,:),winPartialCorLassoSeed(j,:));
            correl3=corrcoef(optimalWinPartialCorLassoSeed(j,:),winPartialCorLassoSeed(j,:));
            correlEachSeedType1(j,i,k)=correl1(1,2);
            correlEachSeedType2(j,i,k)=correl2(1,2);
            correlEachSeedType3(j,i,k)=correl3(1,2);
            
            % use euclidean to estimate distance between FC windows
            % estimated using different lambda
            disp(['Computing euclidean distance for window ',num2str(j)])
            dist1=sqrt(sum(winFullCorSeed(j,:)-winFullCorLassoSeed(j,:)).^2);
            dist2=sqrt(sum(winFullCorSeed(j,:)-winPartialCorLassoSeed(j,:)).^2);
            dist3=sqrt(sum(optimalWinPartialCorLassoSeed(j,:)-winPartialCorLassoSeed(j,:)).^2);
            distEachSeedType1(j,i,k)=dist1;
            distEachSeedType2(j,i,k)=dist2;
            distEachSeedType3(j,i,k)=dist3;
            
            % use mutual information to estimate similarity
            % between FC window creatd using different lambda
            disp(['Computing mutual information for window ',num2str(j)])
            mutualInfo1=mi(winFullCorSeed(j,:),winFullCorLassoSeed(j,:));
            mutualInfo2=mi(winFullCorSeed(j,:),winPartialCorLassoSeed(j,:));
            mutualInfo3=mi(optimalWinPartialCorLassoSeed(j,:),winPartialCorLassoSeed(j,:));
            mutualInfoEachSeedType1(j,i,k)=mutualInfo1;
            mutualInfoEachSeedType2(j,i,k)=mutualInfo2;
            mutualInfoEachSeedType3(j,i,k)=mutualInfo3;
        end
    end
end

similarityMetric={'correl','dist','mutualInfo'};
numMetric=length(similarityMetric);
compareType={'Type1','Type2','Type3'};
numType=length(compareType);
for m=1:numMetric
    metric=char(similarityMetric{m});
    for n=1:numType
        Type=char(compareType{n});
        save([resultDir,metric,'EachSeed',Type,'.mat'], [metric,'EachSeed',Type])
    end
end
end
