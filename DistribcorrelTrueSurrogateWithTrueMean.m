clear
clc
close all


dataLength='all_10min';
session='session1';

resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength, filesep, session, filesep];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/', dataLength,filesep, session,'/surrogateData/'];
tmp1=load([resultDir,'/clustMean/clusterMean_5clusters_session1_normWin.mat']);
clustMean=tmp1.finalMeanWinOfClust;

tmp2=load([resultDir,'zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'.mat']);
winAllSubAllSeedFullCor=tmp2.zWinFullCorLasso;
numROI=156;
numWin=size(winAllSubAllSeedFullCor, 1);
indx=load([resultDir, 'clustIndxNormWinAllSeeds_FullCorLasso_',session,'_10min.txt']);
numClust=length(unique(indx));

disp(['Working on ',session,'_',num2str(numClust),' clusters.'])

corAllStates=zeros(5620, numClust); %5620 is the largest number of windows one state have
for i=1:numClust
    indxList=find(indx==i);
    numWinInClust=length(indxList);
    allWinInClust=zeros(numWinInClust,numROI);
    for j=1:numWinInClust
        for m=1:numROI
            allWinInClust(j,m)=winAllSubAllSeedFullCor(indxList(j), m);
        end
        [cor1, p1]=corrcoef(squeeze(clustMean(i, :)), squeeze(allWinInClust(j, :)));
        corAllStates(j,i)=cor1(1,2);
    end
end

disp('Correlation of real data is done!')

numSurrogate=100;

corAllStatesSurrogate=zeros(5620, numClust,numSurrogate);

for k=1:numSurrogate
    disp(['Working on surrogate ', num2str(k)])
    tmp3=load([resultDir,'surrogate/zNormWindows/zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'_surrogate_', num2str(k),'.mat']);
    winAllSubAllSeedFullCorSurrogate=tmp3.zWinFullCorLasso;
    for i=1:numClust
        indxList=find(indx==i);
        numWinInClust=length(indxList);
        allWinInClustSurrogate=zeros(numWinInClust,numROI);
        for j=1:numWinInClust
            for m=1:numROI
                allWinInClustSurrogate(j,m)=winAllSubAllSeedFullCorSurrogate(indxList(j), m);
            end
            [cor2, p2]=corrcoef(squeeze(clustMean(i, :)), squeeze(allWinInClustSurrogate(j, :)));
            corAllStatesSurrogate(j,i, k)=cor2(1,2);
        end
    end
end

save([resultDir,'surrogate/corWinFCand ClustMeanSurrogate.mat'], 'corAllStatesSurrogate') 

tmp=load([resultDir,'surrogate/corWinFCand ClustMeanSurrogate.mat']);
corAllStatesSurrogate=tmp.corAllStatesSurrogate;

corAllStates1D=reshape(corAllStates, 1, []);
corAllStates1D(find(corAllStates1D==0))=[];
length(corAllStates1D)

corAllStatesSurrogate1D=reshape(corAllStatesSurrogate, 1, []);
corAllStatesSurrogate1D(find(corAllStatesSurrogate1D==0))=[];
length(corAllStatesSurrogate1D)

nElements=hist(corAllStates1D, -1:0.05:1);
norm1=100*nElements/(numWin);

nElements1=hist(corAllStatesSurrogate1D, -1:0.05:1);
norm2=100*nElements1/(numWin*numSurrogate);

close all
figure(1)
bar(norm1, 'r')
hold on
bar(norm2)
set(gca, 'xTick', 1:4:41)
set(gca, 'xTickLabel', -1:0.2:1)
xlabel('Correlation coefficient')
ylabel('Percent ')
saveas(figure(1),[figDir,'session1_corWindFCClustMean_RealSurrogate.png'])



