clear
clc

% compute and plot logliklihood

seed_ROISignals = load('/Users/zhenyang/Documents/microstate/data/time_series/ROISignals_seed_ROISignal.mat');
TC1=seed_ROISignals.ROISignals;
N_seed=size(TC1,2);
%N_seed=1;
ROI_ROISignals=load('/Users/zhenyang/Documents/microstate/data/time_series/ROISignals_atlas_ROISignal.mat');
TC2=ROI_ROISignals.ROISignals;
N_ROI=size(TC2,2);
TC_seed1=[TC1(:,1),TC2];
TC_seed2=[TC1(:,2),TC2];
TC_seed3=[TC1(:,3),TC2];
TC_seed4=[TC1(:,4),TC2];

first_half=TC_seed1(1:size(TC_seed1,1)/2,:);
lambdaList=[5, 10, 20, 50, 100, 200];
lambda=5;
[w, theta, iter, avgTol, hasError] = GraphicalLasso(first_half, lambda)
%[wList, thetaList, lambdaList, errors] = GraphicalLassoPath(first_half, lambdaList)