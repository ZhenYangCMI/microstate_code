clear
clc

% the columns in this file represent sub_ID, DKEFS_DESIGN_SWITCHING_SCALED	DKEFS_VERB_CATEGORYSWITCH_SCALED	DKEFS_TOTALWEIGHTACHIEVESCORE_SCALED	DKEFS_COLORWORD_INHIBITION_SCALED	WASI_PERF	WASI_VERB	WASI_FULL_4
fileDir='/home/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/';
figDir='/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/2sessions/'

a=load([fileDir, 'Neuropsych16sub.txt']);
figure(1)
boxplot(a(:,2:5), 'whisker',3)
saveas(figure(1), [figDir, 'boxplot_DKEFS.png'])
figure(2)
boxplot(a(:,6:8), 'whisker',3)
saveas(figure(1), [figDir, 'boxplot_IQ.png'])


% the metrics are avg_PrctDifNumState	avg_TCPctOverlap	avg_subSeedPrct
% avg_stateFreq	avg_duration	avg_transitions. The sub_ID see sheet
% EmpDyn_session1_2sesCmb21sub.
b=load([fileDir, 'avg_dynamic_metrics.txt']);
numMetric=size(b, 2)
for i=1:numMetric
 figure(i)
boxplot(b(:,i), 'whisker', 3)
saveas(figure(i), [figDir, 'boxplot_avg_dynamic_metric', num2str(i), '.png'])
end