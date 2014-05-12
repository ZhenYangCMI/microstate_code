function [finalWin]=createWin(winWidth, step, TR, numVol)
% output a convolved window
% input
% 1.winWidth: in TRs
% 2.step: sliding steps in TRs
% 3.TR: in seconds
% 4.numVol: totoal time points in TRs

% time
T=0:winWidth-1;

% sampling rate
fs=1/TR;

% number of windows
numWin=floor((numVol-winWidth)/step)+1;

% creat the rectangular window
rectangWin = ones(1,winWidth);

% Calculate gaussian
% Gaussian Mean
GaussMean=(winWidth-1)/2;
% Gaussian Standard Deviation in TRs
GaussStdv=3;
GaussWin=exp(-(T-GaussMean).^2/(2*GaussStdv.^2));

% convolve rectangular and gaussian
convWin=conv(rectangWin,GaussWin);

% length(convWin) = length(rectangWin)+length(GaussWin)-1
% the length(finalWin)=winWidth and take the middle section of the convWin
finalWin=convWin(1,ceil(winWidth/2):(winWidth+ceil(winWidth/2)-1))';

% plot the windows
figure;
subplot(2,2);
subplot(2,2,1); plot(T,rectangWin)
xlim([0 length(rectangWin)-1])
subplot(2,2,2); plot(GaussWin)
xlim([0 length(gaussWin)-1])
subplot(2,2,3); plot(convWin)
xlim([0 length(convWin)-1])
subplot(2,2,4);plot(finalWin)
xlim([0 length(finalWin)-1])