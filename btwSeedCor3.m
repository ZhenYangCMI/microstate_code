clear
clc
close all

numWinPerSeed=5984;
numSeed=4;
figDir='/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/GSR/session1/';
dataDir='/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/GSR/session1/';

% load the iFC windows of all seeds
load([dataDir, 'winFullCorLassoSeed_OptimalLambdaPerSub_645_session1win69_GSR.mat'])
tmp=winFullCorLassoSeed;

% load the index file
indx=load([dataDir, 'clustIndxNormWinAllSeeds_FullCorLasso_session1_10min.txt']);

% recode state 1 and state 2 to match the figs in the paper
indxRecode=indx;
indxRecode(indx==1)=2;
indxRecode(indx==2)=1;

b=[];
preferedState=[2,4,1,1];
for i=1:numSeed
disp(['seed', num2str(i)])
            a=indxRecode(1+numWinPerSeed*(i-1):numWinPerSeed*i);
        for j=1:5984
            if a(j)==preferedState(i)
                a(j)=1;
            else
                a(j)=0;
            end
        end
b(i,:)=a;
end
overlap=sum(b);


seed1=tmp(1+5984*(1-1):5984*1, :);
seed2=tmp(1+5984*(2-1):5984*2, :);
seed3=tmp(1+5984*(3-1):5984*3, :);
seed4=tmp(1+5984*(4-1):5984*4, :);
FC=zeros(length(preferedState),1);

for m=1:length(preferedState)
numWinIncluded=length(find(overlap==m))/5984
seed1cor=seed1(find(overlap==m), 2:4);
seed2cor=[seed2(find(overlap==m),1),seed2(find(overlap==m),3:4)];
seed3cor=[seed3(find(overlap==m),1:2),seed3(find(overlap==m),4)];
seed4cor=[seed4(find(overlap==m),1:3)];

seed1avg=sum(sum(seed1cor))/(length(find(overlap==m))*3);
seed2avg=sum(sum(seed2cor))/(length(find(overlap==m))*3);
seed3avg=sum(sum(seed3cor))/(length(find(overlap==m))*3);
seed4avg=sum(sum(seed4cor))/(length(find(overlap==m))*3);

avgFC(m,1)=seed1avg+seed2avg+seed3avg+seed4avg;
end


