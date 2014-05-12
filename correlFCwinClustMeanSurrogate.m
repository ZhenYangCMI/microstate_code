clear
clc


numWinList=[5602, 5620, 4822,4632, 3260];
numClust=5;
numROI=156;
session='session1'
numSurrogate=100

resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min', filesep, session,filesep];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min',filesep, session,'/surrogateData/'];

clustIndxList=randperm(23936);
    clustIndx1=clustIndxList(1:5602);
    clustIndx2=clustIndxList(5603:5602+5620);
    clustIndx3=clustIndxList(5603+5620+1:5602+5620+4822);
    clustIndx4=clustIndxList(5603+5620+4822+1:5602+5620+4822+4632);
    clustIndx5=clustIndxList(5603+5620+4822+4632+1:5602+5620+4822+4632+3260);
    

% finalMeanWinOfclust=zeros(numClust, numROI);
% for k=1:numSurrogate
%     disp(['Working on surrogate ', num2str(k)])
%         tmp=load([resultDir,'surrogate/zNormWindows/zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'_surrogate_',num2str(k),'.mat']);
%     winAllSubAllSeedFullCor=tmp.zWinFullCorLasso;
%     for i=1:numClust
%         if i==1
%             clustIndx=clustIndx1;
%         elseif i==2
%             clustIndx=clustIndx2;
%         elseif i==3
%             clustIndx=clustIndx3;
%         elseif i==4
%             clustIndx=clustIndx4;
%         else
%             clustindx=clustIndx5;
%         end
%         
%         numWinInClust=length(clustIndx)
%         allWinInClust=zeros(numWinInClust,numROI);
%         for j=1:numWinInClust
%             for m=1:numROI
%                 allWinInClust(j,m)=winAllSubAllSeedFullCor(clustIndx(j), m);
%             end
%         end
%         meanWinOfClust=mean(allWinInClust);
%         for m=1:numROI
%             finalMeanWinOfClust(i,m)=meanWinOfClust(1,m);
%         end
%     end
%     disp('Cluster centroids computed. Start saving file.')
%     save([resultDir, sprintf('surrogate/clustMean/clusterMean_%dclusters_%s_normWin_surrogate_%d.mat',numClust,session,k)],'finalMeanWinOfClust')
% end


corAllStates=zeros(5620, numClust, numSurrogate); %5620 is the largest number of windows one state have
for k=1:numSurrogate
    tmp1=load([resultDir, sprintf('surrogate/clustMean/clusterMean_%dclusters_%s_normWin_surrogate_%d.mat',numClust,session,k)]);
    clustMean=tmp1.finalMeanWinOfClust;
    for i=1:numClust
        
        if i==1
            clustIndx=clustIndx1;
        elseif i==2
            clustIndx=clustIndx2;
        elseif i==3
            clustIndx=clustIndx3;
        elseif i==4
            clustIndx=clustIndx4;
        else
            clustindx=clustIndx5;
        end
        
        numWinInClust=length(clustIndx)
        allWinInClust=zeros(numWinInClust,numROI);
        for j=1:numWinInClust
            for m=1:numROI
                allWinInClust(j,m)=winAllSubAllSeedFullCor(clustIndx(j), m);
            end
            [cor1, p1]=corrcoef(squeeze(clustMean(i, :)), squeeze(allWinInClust(j, :)));
            corAllStates(j,i, k)=cor1(1,2);
        end
    end
end
save([resultDir, sprintf('surrogate/clustMean/corAllStates_%dclusters_%s_normWin.mat',numClust,session)],'corAllStates')
tmp=load([resultDir, sprintf('surrogate/clustMean/corAllStates_%dclusters_%s_normWin.mat',numClust,session)]);
corAllStates=tmp.corAllStates;
corAllStates1DClust1=reshape(corAllStates(1:5602, 1,:),1,[]);
corAllStates1DClust2=reshape(corAllStates(1:5620, 2,:),1,[]);
corAllStates1DClust3=reshape(corAllStates(1:4822, 3,:),1,[]);
corAllStates1DClust4=reshape(corAllStates(1:4632, 4,:),1,[]);
corAllStates1DClust5=reshape(corAllStates(1:3260, 5,:),1,[]);
corAllStates1D=[corAllStates1DClust1, corAllStates1DClust2, corAllStates1DClust3, corAllStates1DClust4, corAllStates1DClust5];
length(corAllStates1D)

tmp1=load([resultDir,'/clustMean/clusterMean_5clusters_session1_normWin.mat']);
clustMean1=tmp1.finalMeanWinOfClust;
tmp2=load([resultDir,'zWinFullCorLasso_OptimalLambdaPerSub_645_',session,'.mat']);
winAllSubAllSeedFullCor1=tmp2.zWinFullCorLasso;
indx=load([resultDir, 'clustIndxNormWinAllSeeds_FullCorLasso_',session,'_10min.txt']);
corAllStatesReal=zeros(5620, numClust); %5620 is the largest number of windows one state have
for i=1:numClust
    indxList=find(indx==i);
    numWinInClust=length(indxList);
    allWinInClust=zeros(numWinInClust,numROI);
    for j=1:numWinInClust
        for m=1:numROI
            allWinInClust(j,m)=winAllSubAllSeedFullCor1(indxList(j), m);
        end
        [cor2, p2]=corrcoef(squeeze(clustMean1(i, :)), squeeze(allWinInClust(j, :)));
        corAllStatesReal(j,i)=cor2(1,2);
    end
end

corAllStatesReal1D=reshape(corAllStatesReal, 1, []);
corAllStatesReal1D(find(corAllStatesReal1D==0))=[];
length(corAllStatesReal1D)

numWin=23936
nElements=hist(corAllStates1D, -1:0.05:1);
norm1=100*nElements/(numWin*numSurrogate);

nElements1=hist(corAllStatesReal1D, -1:0.05:1);
norm2=100*nElements1/numWin;

close all
figure(1)
bar(norm2, 'r')
hold on
bar(norm1)
set(gca, 'xTick', 1:4:41)
set(gca, 'xTickLabel', -1:0.2:1)
xlabel('Correlation coefficient')
ylabel('Percent ')
saveas(figure(1),[figDir,'session1_corWindFCClustMean_RealSurrogate_ownClustMean.png'])



