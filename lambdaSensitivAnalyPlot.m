function [totalNumFigs]=lambdaSensitivAnalyPlot(TR,session, sub)
% This function plot the similarity/distance (estimated using 3 metrics: correlation, euclidean distance, mutual information) between the matrix of
% interests (covariance matrix W and inverse covraince matrix estimated
% suing GLasso) and the base matrix (full correlation estimated with matlab
% correcoef function or the inverse covariance matrix estimated at the
% optimal lambda using GLasso

% clear
% clc
% close all
% 
% sub='8574662';
% session='session1'
% TR={'645'};

lambdaList=(0.08:0.01:0.14);
numLambda=length(lambdaList)
numSeed=4;
numWin=128;
numWinAllSeed=numWin*numSeed;

analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];
resultDir=[analyDir,'results/',char(TR), '/',session,'/lambdaSensitivity/seedFCWin/', sub,'/'];
figDir=[analyDir,'fig/',char(TR), '/',session,'/lambdaSensitivity/', sub,'/'];

%similarityMetric={'correl','dist','mutualInfo'};
similarityMetric={'correl','dist'};
numMetric=length(similarityMetric);
compareType={'Type1','Type2','Type3'};
numType=length(compareType);

% plot the graph for each seed separately
k=0;
for m=1:numMetric
    metric=char(similarityMetric{m});
    for n=1:numType
        Type=char(compareType{n});
        tmp=load([resultDir,metric,'EachSeed',Type,'.mat']);
        data=tmp.([metric,'EachSeed',Type]);
        for i=1:numSeed
            k=k+1
            for j=1:numWin
                dataPlot=data(j,:,i);
                figure(k)
                plot(lambdaList,dataPlot,'Color',[0.8 0.8 0.8],'LineWidth',0.2)
                set(gca,'XTick',0.08:0.01:0.14);
                title([metric,' ', Type,' seed',num2str(i),' ', session,' ',sub]);
                xlabel('Lambda');
                ylabel([metric,' ',Type]);
                xlim([0.08 0.14]);
                if strcmp(metric,'correl') && (strcmp(Type,'Type1') | strcmp(Type,'Type3'))
                    set(gca,'YTick',0.8:0.05:1);
                    ylim([0.8 1]);
                elseif strcmp(metric,'correl') && strcmp(Type, 'Type2')
                    set(gca,'YTick',-0.6:0.05:-0.3);
                    ylim([-0.6 -0.3]);
                elseif strcmp(metric,'dist') && (strcmp(Type,'Type1') | strcmp(Type,'Type2'))
                    set(gca,'YTick',0:5:20);
                    ylim([0 20]);
                else strcmp(metric,'dist') && strcmp(Type, 'Type3')
                    set(gca,'YTick',0:0.2:1);
                    ylim([0 1]);
                end
                hold on
            end
            dataSeed=squeeze(data(:,:,i));
            meanDataPlot=mean(dataSeed);
            plot(lambdaList,meanDataPlot,'b','LineWidth',2)
            hold off
            saveas(figure(k),[figDir,metric,'EachSeed',Type,'Seed',num2str(i),'.png'])
        end
    end
end

% plot the graph for all seeds together
close all
p=0;
for m=1:numMetric
    metric=char(similarityMetric{m});
    for n=1:numType
        Type=char(compareType{n});
        tmp=load([resultDir,metric,'EachSeed',Type,'.mat']);
        data3D=tmp.([metric,'EachSeed',Type]);
        p=p+1
        for i=1:numSeed
            dataSeed=squeeze(data3D(:,:,i));
            meanDataSeed=mean(dataSeed);
            figure (p)
            plot(lambdaList,meanDataSeed,'Color',[0.8 0.8 0.8],'LineWidth',0.6)
            set(gca,'XTick',0.08:0.01:0.14)
            title([metric,' ', Type,' ', session,' ',sub]);
            xlabel('Lambda');
            ylabel([metric,' ',Type]);
            xlim([0.08 0.14])
            if strcmp(metric,'correl') && (strcmp(Type,'Type1') || strcmp(Type,'Type3'))
                set(gca,'YTick',0.8:0.05:1);
                ylim([0.8 1]);
            elseif strcmp(metric,'correl') && strcmp(Type, 'Type2')
                set(gca,'YTick',-0.6:0.05:-0.3)
                ylim([-0.6 -0.3]);
            elseif strcmp(metric,'dist') && (strcmp(Type,'Type1') || strcmp(Type,'Type2'))
                set(gca,'YTick',0:5:20);
                ylim([0 20]);
            else strcmp(metric,'dist') && strcmp(Type, 'Type3')
                set(gca,'YTick',0:0.2:1);
                ylim([0 1]);
            end
            hold on
        end
        data2D=vertcat(data3D(:,:,1),data3D(:,:,2),data3D(:,:,3),data3D(:,:,4));
        meanData=mean(data2D);
        plot(lambdaList,meanData,'b','LineWidth',2)
        hold off
        saveas(figure(p),[figDir,metric,'AllSeeds',Type,'.png'])
    end
end

totalNumFigs=k+p;
end