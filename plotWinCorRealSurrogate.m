clear
clc
close all

numSurrogate=100;
dataLength='all_10min';
session='session1';

resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/GSR/session1/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/', dataLength,'/GSR/', session,'/surrogateData/'];

tmp=load([resultDir, 'winFullCorLasso_OptimalLambdaPerSub_645_session1.mat']);
real=tmp.winFullCorLasso;
real1D=reshape(real,1,[]);
numWin=size(real, 1);
numROI=size(real,2);

surrogateAll=zeros(numWin, numROI, numSurrogate);
for i=1:numSurrogate
    disp(['Working on surrogate ', num2str(i)])
    tmp1=load([resultDir, 'surrogate/unNormWindows/winFullCorLasso_OptimalLambdaPerSub_645_session1_surrogate_', num2str(i),'.mat']);
    surrogate=tmp1.winFullCorLasso;
    surrogateAll(:,:,i)=surrogate;
end

surrogate1D=reshape(surrogateAll, 1,[]);

nElements=hist(real1D, -1:0.05:1);
norm1=100*nElements/(numWin*numROI);

nElements1=hist(surrogate1D, -1:0.05:1);
norm2=100*nElements1/(numWin*numROI*numSurrogate);


% plot for the paper
close all
figure(1)
h1=bar(norm1, 'r', 'LineWidth', 0.5, 'BarWidth', 1)
hold on
h2=bar(norm2, 'FaceColor', [0/255, 0/255,1], 'LineWidth', 0.5, 'BarWidth', 1)
set(gca, 'xTick', 1:4:41)
set(gca, 'xTickLabel', -1:0.2:1)
box off
set(gca, 'LineWidth', 2)
ch = get(h2,'child');
set (ch, 'FaceAlpha', 0.5);
% xlabel('Correlation coefficient')
% ylabel('Percent ')
saveas(figure(1),[figDir,'session1_corOfAllWins_RealSurrogate_paper.png'])

overallWidth=10;
h1=BarSpecial(norm1, overallWidth)

close all
figure(2)
bar(norm1, 'r')
hold on
bar(norm2,'b')
set(gca, 'xTick', 1:4:41)
set(gca, 'xTickLabel', -1:0.2:1)
xlabel('Correlation coefficient')
ylabel('Percent ')
saveas(figure(2),[figDir,'session1_corOfAllWins_RealSurrogate.png'])




