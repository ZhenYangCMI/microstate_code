function [correlAllSeedsType3, distAllSeedsType3,mutualInfoAllSeedsType3]=lambdaSensitivAnalySeedTogether(TR,session, sub, optimalLambda)
%the following function estimate the sensitivity of lambda to the FC
%windows using 3 metric: correlation, euclidean distance, and mutual
%information. There are three types of comparisons across 7 lambdas. Type1: compare the similarity between the full correlation matrix generated using Lasso
%(i.e. covariance matrix computed on normalized TS of that window; Type2:
%compare the full correlation matrix generated using matlab corrcoef function and
%the inverse covariance matrix generated using GLasso on normalized TS of
%that window (i.e. partial correlation); Type 3; compare the parital
%correlation matrix at differen lambda with the partial correlation matrix
%computed at the lambdaAtMaxLog
         

%1 TR: cell, e.g. TR={'645'}; 2mutualInfoAllSeedsType1 session: string, e.g 'session1'; 3.
% sub:string, e.g. '8574662'; 4. optimalLambda: numerical

% clear
% clc
% 
% sub='8574662';
% session='session1'
% TR={'645'};
% optimalLambda=0.12;


lambdaList=(0.08:0.01:0.14);
numLambda=length(lambdaList)
numWin=128;
numSeed=4;
numWinAllSeed=numWin*numSeed
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
    
            for j=1:numWinAllSeed
            % use correlation to estimate similarity
            disp(['Computing correlations for window ',num2str(j)])
            correl1=corrcoef(winFullCor(j,:),winFullCorLasso(j,:));
            correl2=corrcoef(winFullCor(j,:),winPartialCorLasso(j,:));
            correl3=corrcoef(optimalWinPartialCorLasso(j,:),winPartialCorLasso(j,:));
            correlAllSeedsType1(j,i)=correl1(1,2);
            correlAllSeedsType2(j,i)=correl2(1,2);
            correlAllSeedsType3(j,i)=correl3(1,2);
            
            % use euclidean to estimate distance between FC windows
            % estimated using different lambda
            disp(['Computing euclidean distance for window ',num2str(j)])
            dist1=sqrt(sum(winFullCor(j,:)-winFullCorLasso(j,:)).^2);
            dist2=sqrt(sum(winFullCor(j,:)-winPartialCorLasso(j,:)).^2);
            dist3=sqrt(sum(optimalWinPartialCorLasso(j,:)-winPartialCorLasso(j,:)).^2);
            distAllSeedsType1(j,i)=dist1;
            distAllSeedsType2(j,i)=dist2;
            distAllSeedsType3(j,i)=dist3;
            
            % use mutual information to estimate similarity
            % between FC window creatd using different lambda
            disp(['Computing mutual information for window ',num2str(j)])
            mutualInfo1=mi(winFullCor(j,:),winFullCorLasso(j,:));
            mutualInfo2=mi(winFullCor(j,:),winPartialCorLasso(j,:));
            mutualInfo3=mi(optimalWinPartialCorLasso(j,:),winPartialCorLasso(j,:));
            mutualInfoAllSeedsType1(j,i)=mutualInfo1;
            mutualInfoAllSeedsType2(j,i)=mutualInfo2;
            mutualInfoAllSeedsType3(j,i)=mutualInfo3;
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
        save([resultDir,metric,'AllSeeds',Type,'.mat'], [metric,'AllSeeds',Type])
    end
end
end
