


clear
clc
close all


numSeed=4;
numSub=21;
numState=5;


dataLength='all_10min';
load('MyColormaps','mycmap') % [0 1]
load('MyColormaps2','mycmap2') % > 1
load('MyColormaps3','mycmap3') % log scale, white largest value

resultDir='/Users/zhenyang/Desktop/Zhen/results/all_10min/2sesConcate/';
figDir='/Users/zhenyang/Desktop/Zhen/figs/5_3_13/';

measure='probTM'
tmp1=load([resultDir, measure,'_session1.mat']);
data1=tmp1.(measure);
tmp2=load([resultDir, measure,'_session2.mat']);
data2=tmp2.(measure);

measure1='meanProbTM';
tmp3=load([resultDir, measure1,'_session1.mat']);
data3=tmp3.(measure1);
tmp4=load([resultDir, measure1,'_session2.mat']);
data4=tmp4.(measure1);

measure2='meanProbTMDifState';
tmp5=load([resultDir, measure2,'_session1.mat']);
data5=tmp5.(measure2);
tmp6=load([resultDir, measure2,'_session2.mat']);
data6=tmp6.(measure2);

% data1(numStateT, numStateT+1, numSub, numSeed)

%Seed 2, stateT=1, stateT+1=1
seed=2;
stateT=2;
stateT1=4;
a=squeeze(data1(stateT,stateT1,:,seed));
b=squeeze(data2(stateT,stateT1,:,seed));
sesCmb=[a,b];
ses1NumZeros=length(find(a==0))
ses2NumZeros=length(find(b==0))
%save(sprintf('/Users/zhenyang/Desktop/Zhen/results/all_10min/2sesConcate/TM_seed%dState%dToState%d.mat',seed,stateT,stateT1), 'sesCmb')