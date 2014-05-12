function [finalWin]=winCreation(winWidth, plotWin)
% output the middle part of a convolved window as the final win applied to
% TS
% input:
% 1.winWidth: in TRs
% 2.plotWin: plot the win or not (1=yest, 0=no)


% time
T=0:winWidth-1;

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
finalWin=convWin(1,ceil(winWidth/2):(winWidth+ceil(winWidth/2)-1));

if (plotWin==1)
    % plot the windows
    figure;
    subplot(2,2,1); plot(T,rectangWin)
    xlim([0 length(rectangWin)-1])
    subplot(2,2,2); plot(GaussWin)
    xlim([0 length(GaussWin)-1])
    subplot(2,2,3); plot(convWin)
    xlim([0 length(convWin)-1])
    subplot(2,2,4);plot(finalWin)
    xlim([0 length(finalWin)-1])
else
    disp('Final window created! It will not be plotted.')
end