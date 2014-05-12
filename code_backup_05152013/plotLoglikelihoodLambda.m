function [ output_args ] = logLambdaPlot( lambdaList, loglikelihoodList, sub)
% plot loglikelihood as a function of lambda and save the figs
%

figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/'];
for q=1:numWin
    figure(1)
    subplot(2,3,q)
    plot(lambdaList(q,:), logLikelihoodList(q,:))
end

saveas(figure(1),[figDir,'lamdaEstimate_',session,'_sub_', char(subID), '.png'])

end

