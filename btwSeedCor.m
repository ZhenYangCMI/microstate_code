clear
clc

% load the iFC windows and define x as windows of seed 1 and y as windows
% of seed 2
load('winFullCorLassoSeed_OptimalLambdaPerSub_645_session1win69_GSR.mat')
tmp=winFullCorLassoSeed;
x=tmp(1:5984, :);
y=tmp(5985:5984*2, :);

% load the index file and assign index of seed 1 to indxX and index of seed
% 2 to indxY
indx=load('clustIndxNormWinAllSeeds_FullCorLasso_session1_10min.txt');
indxX=indx(1:5984);
indxY=indx(5985:5984*2);

% binarize the index and assign it to xTime and yTime
xTime(indxX==1)=1;
yTime(indxY==4)=1;
yTime(5901:5984)=0; % The last window with indxY==4 is 5900, so zeros were added to make yTime the same length

% z is vector with 0 (two seeds are not in oppostive states) or 1 (two seeds are in the opposite states)
z=(yTime == xTime).* yTime;

opp=find(z);
numOpp=length(opp)
notOpp=find(z==0);
numNotOpp=length(notOpp)

%% when the two seeds are in oppositve states

% The average correlations between seed1/seed2 and other seeds (seed 3 and
% seed 4 with themselves excluded)
a=x(opp, 3)+x(opp,4)+y(opp,3)+y(opp,4);
aAvg=sum(a)/(numOpp*4)

% The average correlations between seed 1 and other seeds (seed 3 and seed 4)
a1=x(opp, 3)+x(opp,4);
a1Avg=sum(a1)/(numOpp*2)

% The average correlations between seed 2 and other seeds (seed 3 and seed 4)
a2=y(opp,3)+y(opp,4);
a2Avg=sum(a2)/(numOpp*2)

% The average correlations between seed 1 and seed 2 
b=x(opp, 2)+y(opp,1);
bAvg=sum(b)/(numOpp*2)

%% when the two seeds are not in oppositve states

% The average correlations between seed1/seed2 and other seeds (seed 3 and
% seed 4 with themselves excluded)
c=x(notOpp, 3)+x(notOpp,4)+y(notOpp,3)+y(notOpp,4);
cAvg=sum(c)/(numNotOpp*4)

% The average correlations between seed 1 and other seeds (seed 3 and seed 4)
c1=x(notOpp, 3)+x(notOpp,4);
c1Avg=sum(c1)/(numNotOpp*2)

% The average correlations between seed 2 and other seeds (seed 3 and seed 4)
c2=y(notOpp,3)+y(notOpp,4);
c2Avg=sum(c2)/(numNotOpp*2)

% The average correlations between seed 1 and seed 2 
d=x(notOpp, 2)+y(notOpp,1);
dAvg=sum(d)/(numNotOpp*2)


