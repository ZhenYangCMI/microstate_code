
clear
clc
close all

numSeed=4;
numSub=22;
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/session1/ses1_ward_euclidean/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/session1/ses1_ward_euclidean/'];

for i=2:2:16
    numClust=i;
    tmp=load([resultDir,'indx_',num2str(numClust),'clusters.mat']);
    indx=tmp.indx;
    [seedPrct, subPrct, subSeedPrct] = DoQuickStats(indx,numClust);
    save([resultDir, 'seedPrct_',num2str(numClust),'clusters.mat'],'seedPrct')
    save([resultDir, 'subPrct_',num2str(numClust),'clusters.mat'],'subPrct')
    save([resultDir, 'subSeedPrct_',num2str(numClust),'clusters.mat'],'subSeedPrct')
end

%categoryList={'seedPrct', 'subPrct','subSeedPrct'};
categoryList={'subSeedPrct'};
numCategory=length(categoryList);
for j=1:numCategory
    k=0;
    categoryName=char(categoryList{j});
    for i=2:2:16
        numClust=i;
        disp( ['Working on ',categoryName, '_',num2str(numClust)])
        k=k+1;
        tmp=load([resultDir,categoryName,'_',num2str(numClust),'clusters.mat']);
        category=tmp.(categoryName);
        clear dataPlot1 dataPlot2
        if strcmp('subSeedPrct',categoryName)
            for m=1:numSub
                dataPlot1=squeeze(category(:,:,m));
                if (m<5)
                    figure(1)
                    subplot(2,2,m)
                    imagesc(dataPlot1)
                    colorbar
                    xlabel('Seeds')
                    ylabel('clusters')
                elseif (m<9)
                    figure(2)
                    subplot(2,2,m-4)
                    imagesc(dataPlot1)
                    colorbar
                    xlabel('Seeds')
                    ylabel('clusters')
                elseif (m<13)
                    figure(3)
                    subplot(2,2,m-8)
                    imagesc(dataPlot1)
                    colorbar
                    xlabel('Seeds')
                    ylabel('clusters')
                elseif (m<17)
                    figure(4)
                    subplot(2,2,m-12)
                    imagesc(dataPlot1)
                    colorbar
                    xlabel('Seeds')
                    ylabel('clusters')
                elseif (m<21)
                    figure(5)
                    subplot(2,2,m-16)
                    imagesc(dataPlot1)
                    colorbar
                    xlabel('Seeds')
                    ylabel('clusters')
                else
                    figure(6)
                    subplot(2,2,m-20)
                    imagesc(dataPlot1)
                    colorbar
                    xlabel('Seeds')
                    ylabel('clusters')
                end
            end
            saveas(figure(1),[figDir,'subSeedPrct/',categoryName,'_plotPerSeed_1-4Sub_',num2str(numClust),'clusters.png'])
            saveas(figure(2),[figDir,'subSeedPrct/',categoryName,'_plotPerSeed_5-8Sub_',num2str(numClust),'clusters.png'])
            saveas(figure(3),[figDir,'subSeedPrct/',categoryName,'_plotPerSeed_9-12sub_',num2str(numClust),'clusters.png'])
            saveas(figure(4),[figDir,'subSeedPrct/',categoryName,'_plotPerSeed_13-16Sub_',num2str(numClust),'clusters.png'])
            saveas(figure(5),[figDir,'subSeedPrct/',categoryName,'_plotPerSeed_17-20Sub_',num2str(numClust),'clusters.png'])
            saveas(figure(6),[figDir,'subSeedPrct/',categoryName,'_plotPerSeed_21-22sub_',num2str(numClust),'clusters.png'])
            
            for n=1:numSeed
                dataPlot2=squeeze(category(:,n,:));
                
                if (n<3)
                    figure(7)
                    subplot(2,1,n)
                    imagesc(dataPlot2)
                    colorbar
                    title('Percent of windows per subject')
                    xlabel('Subjects')
                    ylabel('clusters')
                else
                    figure(8)
                    subplot(2,1,n-2)
                    imagesc(dataPlot2)
                    colorbar
                    title('Percent of windows per subject')
                    xlabel('Subjects')
                    ylabel('clusters')
                end
                saveas(figure(7),[figDir,'subSeedPrct/',categoryName,'_plotPerSub_first2Seeds_',num2str(numClust),'clusters.png'])
                saveas(figure(8),[figDir,'subSeedPrct/',categoryName,'_plotPerSub_last2Seeds_',num2str(numClust),'clusters.png'])
            end
        else
            figure(k)
            imagesc(category)
            colorbar
            title('Percent of windows')
            xlabel('Seeds')
            ylabel('clusters')
            saveas(figure(k),[figDir,categoryName,'_',num2str(numClust),'clusters.png'])
        end
    end
end
