clear
clc

dataDir=['/home/data/Projects/microstate/DPARSF_preprocessed/data/645/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/'];
motionData=load([dataDir,'motion.txt']);
concateMotion=[motionData(1,:),motionData(2,:)];

figure(1)
subplot(1,2,1)
plot(concateMotion)
ylim([0.06 0.4])
ylabel('Mean FD')
xlabel('Session 1 and Session 2 of the Same Group of Subjects')
subplot(1,2,2)
hist(concateMotion,0.05:0.05:0.5)
ylabel('Frequency')
xlabel('Mean FD')
saveas(figure(1),[figDir, 'motionDistribution.png']);

for i=1:2
    figure(i+1)
    subplot(1,2,1)
    plot(motionData(i,:))
    ylabel('Mean FD')
    xlabel('Subjects')
    ylim([0.06 0.4])
    subplot(1,2,2)
    hist(motionData(i,:),0:0.05:0.5)
    axis([0 0.5 0 12])
    ylabel('Frequency')
    xlabel('Mean FD')
    saveas(figure(i+1),[figDir, 'motionDistributionSes',num2str(i),'.png']);
end

lambda=load('/home/data/Projects/microstate/DPARSF_preprocessed/results/645/meanLambdaAtMaxLogAcrossWin.mat');
lambdaAtMaxLog=lambda.meanLambdaAtMaxLogAcrossWin;
concateLambdaAtMaxLog=[lambdaAtMaxLog(1,:),lambdaAtMaxLog(2,:)];
figure(4)
subplot(1,2,1)
plot(concateLambdaAtMaxLog)
ylim([0.07 0.15])
ylabel('Mean lambdaAtMaxLog')
xlabel('Session 1 and Session 2 of the Same Group of Subjects')
subplot(1,2,2)
hist(concateLambdaAtMaxLog,0.05:0.01:0.18)
axis([0.05 0.18 0 14])
ylabel('Frequency')
xlabel('Mean ambdaAtMaxLog')
saveas(figure(4),[figDir, 'meanLambdDistribution.png']);


figure(5)
subplot(1,2,1)
plot(lambdaAtMaxLog(1,:))
ylim([0.07 0.15])
ylabel('Mean lambdaAtMaxLog')
xlabel('Subjects Session1')
subplot(1,2,2)
hist(lambdaAtMaxLog(1,:),0.05:0.01:0.18)
axis([0.05 0.18 0 14])
ylabel('Frequency')
xlabel('Mean lambdaAtMaxLog')
saveas(figure(5),[figDir, 'lambdaAtMaxLogDistributionSes1.png']);

figure(6)
subplot(1,2,1)
plot(lambdaAtMaxLog(2,:))
ylim([0.07 0.15])
ylabel('Mean lambdaAtMaxLog')
xlabel('Subjects Session2')
subplot(1,2,2)
hist(lambdaAtMaxLog(2,:),0.05:0.01:0.18)
axis([0.05 0.18 0 14])
ylabel('Frequency')
xlabel('Mean lambdaAtMaxLog')
saveas(figure(6),[figDir, 'lambdaAtMaxLogDistributionSes2.png']);



[cor2Sessions, p2Session]=corrcoef(concateMotion,concateLambdaAtMaxLog);
[corSession1,pSession1]=corrcoef(motionData(1,:), lambdaAtMaxLog(1,:));
[corSession2,pSession2] =corrcoef(motionData(2,:), lambdaAtMaxLog(2,:));
[corMotion,pMotion]=corrcoef(motionData(1,:), motionData(2,:));
[corLambda,pLambda] =corrcoef(lambdaAtMaxLog(1,:), lambdaAtMaxLog(2,:));
figure(7)
scatter(motionData(1,:), lambdaAtMaxLog(1,:), 'r')
hold on
scatter(motionData(2,:), lambdaAtMaxLog(2,:),'b')
xlabel('Motion meanFD')
ylabel('LambdaAtMaxLog')
axis([0.05 0.4 0.06 0.16])
lsline

saveas(figure(7),[figDir, 'CorrelMeanFD&Lambda.png']);

figure(8)
subplot(1,2,1)
scatter(motionData(1,:), motionData(2,:))
xlabel('Session 1 Motion meanFD')
ylabel('Session 2 Motion meanFD')
lsline
subplot(1,2,2)
scatter(lambdaAtMaxLog(1,:), lambdaAtMaxLog(2,:))
xlabel('Session 1 LambdaAtMaxLog')
ylabel('Session 2 LambdaAtMaxLog')
lsline
saveas(figure(8),[figDir, 'CorrelMeanFD&LambdaBtwSessions.png']);