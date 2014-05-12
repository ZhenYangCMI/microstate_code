clear
clc
close all

sessionList={'session1','session2'};
numSession=length(sessionList);
resultDir1=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/first_5min/'];
resultDir2=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/'];

lambdaAtMaxLog1=load([resultDir1,'meanLambdaAtMaxLogAcrossWin.txt']);
lambdaAtMaxLog2=load([resultDir2,'meanLambdaAtMaxLogAcrossWin.txt']);

% compute correlation and scatterplot the lambdaAtMaxLog across sessions
[corSession1, pSession2]=corrcoef(squeeze(lambdaAtMaxLog1(:,1)),squeeze(lambdaAtMaxLog2(:,1)))
[corSession2, pSession2]=corrcoef(squeeze(lambdaAtMaxLog1(:,2)),squeeze(lambdaAtMaxLog2(:,2)))

for i=1:numSession
    session=char(sessionList{i});
    figure(i)
    scatter(squeeze(lambdaAtMaxLog1(:,i)),squeeze(lambdaAtMaxLog2(:,i)))
    xlabel(['5min ',session])
    ylabel(['10min ',session])
    set(gca,'XTick',0.05:0.05:0.20,'YTick', 0.05:0.05:0.20);
    xlim ([0.05 0.20])
    ylim ([0.05 0.20])
    title(['lambdaAtMaxLog 5minVS10min ', session])
    lsline
    saveas(figure(i),[figDir,'scatterPlot_5minVS10min_lambdaAtMaxLog_',session,'.png'])
end

% plot distribution of 10 min lambdaAtMaxLog averaged across windows

lambdaAtMaxLog10min=vertcat(lambdaAtMaxLog2(:,1),lambdaAtMaxLog2(:,2));
[h, p]=ttest(lambdaAtMaxLog2(:,1),lambdaAtMaxLog2(:,2))

figure(3)
subplot(1,2,1)
plot(lambdaAtMaxLog10min)
ylim([0.07 0.21])
xlim([0 50])
ylabel('LambdaAtMaxLog')
xlabel('10 min both sessions all subs across wins')
set(gca,'XTick',0:10:50,'YTick', 0.07:0.02:0.21);
subplot(1,2,2)
hist(lambdaAtMaxLog10min,0.07:0.01:0.20)
axis([0.07 0.2 0 15])
ylabel('Frequency')
xlabel('Mean ambdaAtMaxLog')
saveas(figure(3),[figDir, 'all_10min/meanLambdDistributionAcrossWin.png']);


% plot distribution of 10 min lambdaAtMaxLog for all windows
tmp=load([resultDir2,'lambdaAtMaxLog_2sessions_22subjects.mat']);
lambdaAtLogWinSep=tmp.lambdaAtMaxLog;
tmp1=lambdaAtLogWinSep(:,:,1)';
tmp2=lambdaAtLogWinSep(:,:,2)';
lambdaAtLogWinSep1d=vertcat(reshape(tmp1,[],1),reshape(tmp2,[],1));

figure(4)
plot(lambdaAtLogWinSep1d)
ylim([0.07 0.21])
set(gca,'XTick',0:100:600,'YTick', 0.07:0.02:0.21);
ylabel('LambdaAtMaxLog')
xlabel('10 min both sessions all subs all wins')
saveas(figure(4),[figDir, 'all_10min/LambdDistributionWinSep.png']);


%% replot distribution of 5 min lambdaAtMaxLog averaged across windows

lambdaAtMaxLog5min=vertcat(lambdaAtMaxLog1(:,1),lambdaAtMaxLog1(:,2));
[h, p]=ttest(lambdaAtMaxLog1(:,1),lambdaAtMaxLog1(:,2))

figure(5)
subplot(1,2,1)
plot(lambdaAtMaxLog5min)
ylim([0.07 0.21])
xlim([0 50])
ylabel('LambdaAtMaxLog')
xlabel('5 min both sessions all subs across wins')
set(gca,'XTick',0:10:50,'YTick', 0.07:0.02:0.21);
subplot(1,2,2)
hist(lambdaAtMaxLog5min,0.07:0.01:0.20)
axis([0.07 0.2 0 15])
ylabel('Frequency')
xlabel('Mean ambdaAtMaxLog')
saveas(figure(5),[figDir, 'first_5min/meanLambdDistributionAcrossWin.png']);


% replot distribution of 5 min lambdaAtMaxLog for all windows
tmp=load([resultDir1,'lamdaAtMaxLog.mat']);
lambdaAtLogWinSep=tmp.lambdaAtMaxLog;
tmp1=lambdaAtLogWinSep(:,:,1)';
tmp2=lambdaAtLogWinSep(:,:,2)';
lambdaAtLogWinSep1d=vertcat(reshape(tmp1,[],1),reshape(tmp2,[],1));

figure(6)
plot(lambdaAtLogWinSep1d)
ylim([0.07 0.21])
xlim([0 600])
set(gca,'XTick',0:100:600,'YTick', 0.07:0.02:0.21);
ylabel('LambdaAtMaxLog')
xlabel('5 min both sessions all subs all wins')
saveas(figure(6),[figDir, 'first_5min/LambdDistributionWinSep.png']);