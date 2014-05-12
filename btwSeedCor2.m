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

colorList=['b','r','g'];

y=zeros(numWinPerSeed, 1);
for i=1:numSeed
figure(i)
stateIndx=indxRecode(1+numWinPerSeed*(i-1):numWinPerSeed*i);
plot(stateIndx, '-*')
ylim([-1 6])
t=0;
for j=1:numSeed
if j~=i
t=t+1;
color=colorList(t);
seedCor=tmp(1+numWinPerSeed*(i-1):numWinPerSeed*i, j);
hold on
plot(seedCor, color)
end
end
hold on
plot(y, '--k')
saveas(figure(i), [figDir, 'stateTimeCoursAndCorBtwSeeds_seed', num2str(i), '.png'])
end

for i=1:numSeed
seed=[1,2,3,4];
seedAvg=seed;
seedAvg(seed==i)=[];
seedCor=tmp(1+numWinPerSeed*(i-1):numWinPerSeed*i, seedAvg);
seedAvgCor=mean(seedCor, 2);
allSeedAvgCor(1+numWinPerSeed*(i-1):numWinPerSeed*i,1)=seedAvgCor;
end

for i=1:numSeed
   figure(i+4)
stateIndx=indxRecode(1+numWinPerSeed*(i-1):numWinPerSeed*i);
plot(stateIndx, '-*')
ylim([-1 6])
avgCor=allSeedAvgCor(1+numWinPerSeed*(i-1):numWinPerSeed*i);
hold on
plot(avgCor, 'r')
hold on
plot(y, '--k')
saveas(figure(i+4), [figDir, 'stateTimeCoursAndAvgCorAllSeeds_seed', num2str(i), '.png'])
end


% binarize the index and compute the prc overlap

seed2=indxRecode(1+numWinPerSeed:numWinPerSeed*2);
for j=1:5984
    if seed2(j)==4
        seed2(j)=1;
    else
        seed2(j)=0;
    end
end

t=0;
for i=1:numSeed
disp(['seed', num2str(i)])
    if i~=2
        t=t+1;
        a=indxRecode(1+numWinPerSeed*(i-1):numWinPerSeed*i);
        for j=1:5984
            if a(j)==1 || a(j)==2
                a(j)=1;
            else
                a(j)=0;
            end
        end
        overlap=sum(a.*seed2)
        overlapFinal(t,1)=overlap/5984*100
        z=a.*seed2;
        opp=find(z==1);
        numOpp=length(opp)
        notOpp=find(z==0);
        numNotOpp=length(notOpp)
        
        if i==1
            y=tmp(1+numWinPerSeed:numWinPerSeed*2, :);
            a2=y(opp,3)+y(opp,4);
            a2Avg=sum(a2)/(numOpp*2)
            
            b2=y(notOpp,3)+y(notOpp,4);
            b2Avg=sum(b2)/(numNotOpp*2)
            
            x=tmp(1+numWinPerSeed*(i-1):numWinPerSeed*i, :);
            a=x(opp,3)+x(opp,4);
            aAvg=sum(a)/(numOpp*2)
            
            b=x(notOpp,3)+x(notOpp,4);
            bAvg=sum(b)/(numNotOpp*2)
            
        elseif i==3
            y=tmp(1+numWinPerSeed:numWinPerSeed*2, :);
            a2=y(opp,1)+y(opp,4);
            a2Avg=sum(a2)/(numOpp*2)
            
            b2=y(notOpp,1)+y(notOpp,4);
            b2Avg=sum(b2)/(numNotOpp*2)
            
            x=tmp(1+numWinPerSeed*(i-1):numWinPerSeed*i, :);
            a=x(opp,1)+x(opp,4);
            aAvg=sum(a)/(numOpp*2)
            
            b=x(notOpp,1)+x(notOpp,4);
            bAvg=sum(b)/(numNotOpp*2)
            
        else
            y=tmp(1+numWinPerSeed:numWinPerSeed*2, :);
            a2=y(opp,3)+y(opp,1);
            a2Avg=sum(a2)/(numOpp*2)
            
            b2=y(notOpp,3)+y(notOpp,1);
            b2Avg=sum(b2)/(numNotOpp*2)
            
            x=tmp(1+numWinPerSeed*(i-1):numWinPerSeed*i, :);
            a=x(opp,3)+x(opp,1);
            aAvg=sum(a)/(numOpp*2)
            
            b=x(notOpp,3)+x(notOpp,1);
            bAvg=sum(b)/(numNotOpp*2)
        end
    end
end

