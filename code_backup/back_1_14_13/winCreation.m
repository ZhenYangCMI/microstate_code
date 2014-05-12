clear
clc
close all

% window width in TRs
winWidth = 69;
% window sliding step in TRs
step=3;
% time
T=0:winWidth-1;
% repetition time in secs
TR=0.645;
% sampling rate
fs=1/TR;
% total number of time points
numVol=450;
% number of windows
numWin=floor((numVol-winWidth)/step)+1;


% creat the rectangular window
rectangWin = ones(1,winWidth);
% figure(1); subplot(1,2,1); plot(T,RW)
% xlim([0 length(RW)-1])

% Calculate gaussian
% Gaussian Mean
GaussMean=(winWidth-1)/2;
% Gaussian Standard Deviation in TRs
GaussStdv=3;
GaussWin=exp(-(T-GaussMean).^2/(2*GaussStdv.^2));
% figure(2); subplot(1,2,1); plot(GW)
% xlim([0 length(GW)-1])

% convolve rectangular and gaussian
convWin=conv(rectangWin,GaussWin);
% figure(3); subplot(1,2,1); plot(WC)
% xlim([0 length(WC)-1])