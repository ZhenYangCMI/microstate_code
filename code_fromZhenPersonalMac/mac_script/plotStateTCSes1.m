

clear
clc
close all

session='session1'

W_width = 69;
% window sliding step in TRs
step=3;
% total number of time points
N_vol=884;
numSub=22;
numSeed=4;
% number of windows
numWinPerSub=floor((N_vol-W_width)/step)+1;
numWinPerSeed=numWinPerSub*numSub
Time=ceil(W_width/2):step:(ceil(W_width/2)+step*(numWinPerSub-1));

numROI=156;
resultDir=['/Users/zhenyang/Desktop/Zhen/results/all_10min/session1/'];
figDir=['/Users/zhenyang/Desktop/Zhen/figs/6_2_13/'];

indx=load([resultDir,'clustIndxNormWinAllSeeds_FullCorLasso_',session,'_10min.txt']);
numClust=length(unique(indx));
disp ('Files loaded successfully.')

indxRecode=zeros(length(indx),1);

% match the states between two analysis: each session clustered separately and two sessions concatenated
if strcmp(session,'session1')
    indxRecode(find(indx==1))=2;
    indxRecode(find(indx==2))=1;
    indxRecode(find(indx==3))=3;
    indxRecode(find(indx==4))=4;
    indxRecode(find(indx==5))=5;
end

% specify which sub to plot (1, 15, 21)
sub=21
lineWidth1=3
lineWidth2=6
plotcolor=[230,159,0;0,158,115;204,121,167;86,180,233];
plotcolor=plotcolor/255;

for j=sub
    figure(1)
    
    for i=1:numSeed
        
        indxRecodeSeed=indxRecode((1+numWinPerSeed*(i-1)):numWinPerSeed*i);
        indxRecodeSeedSub=indxRecodeSeed((1+numWinPerSub*(j-1)):numWinPerSub*j);
        h=plot(Time,indxRecodeSeedSub,'Color',[plotcolor(i,:)])
        ylim([0 6])
        xlim([0 N_vol])
        set(h,'LineWidth',lineWidth2)
        set(gca,'xTick', 0:200:800,'yTick', 0:6, 'XTickLabel',[],'YTickLabel',[],'box','off')
        set(gca,'LineWidth',lineWidth1)
        %         ylabel('States')
        %         xlabel('Time')
        %         title('Time Course of States')
        %legend('Seed1','Seed2','Seed3','Seed4')
        hold on
    end
    hold off
    saveas(figure(1),[figDir, 'stateTC_sub', num2str(j), '_',session, 'seedsTogether.png'])
end
close all


for j=sub
    for i=1:numSeed
        figure(i)
        indxRecodeSeed=indxRecode((1+numWinPerSeed*(i-1)):numWinPerSeed*i);
        indxRecodeSeedSub=indxRecodeSeed((1+numWinPerSub*(j-1)):numWinPerSub*j);
        h=plot(Time,indxRecodeSeedSub, 'Color',[plotcolor(i,:)])
        set(h,'LineWidth',lineWidth2)
        set(gca,'xTick', 0:200:800,'yTick', 0:6, 'XTickLabel',[],'YTickLabel',[],'box','off')
       set(gca,'LineWidth',lineWidth1)
        ylim([0 6])
        xlim([0 N_vol])
        %         ylabel('States')
        %         xlabel('Time')
        %         title(['Time Course of States -- Seed ', num2str(i)])
        
        saveas(figure(i),[figDir, 'stateTC_sub', num2str(j), '_seed',num2str(i),'_',session, '_seedSep.png'])
    end
end
close all

possibNumStatesAtATime=4; % 4 possibilities
numUniqStat4Seeds=zeros(numWinPerSub, numSub);
prctDifNumStates=zeros(numSub, possibNumStatesAtATime);
for j=sub
    indxSeed1=indxRecode((1+numWinPerSub*(j-1)):numWinPerSub*j);
    (1+numWinPerSub*(j-1))
    numWinPerSub*j
    indxSeed2=indxRecode((1+numWinPerSub*(j-1)+numWinPerSeed):numWinPerSub*j+numWinPerSeed);
    (1+numWinPerSub*(j-1)+numWinPerSeed)
    numWinPerSub*j+numWinPerSeed
    indxSeed3=indxRecode((1+numWinPerSub*(j-1)+numWinPerSeed*2):numWinPerSub*j+numWinPerSeed*2);
    indxSeed4=indxRecode((1+numWinPerSub*(j-1)+numWinPerSeed*3):numWinPerSub*j+numWinPerSeed*3);
    indxAllSeeds=zeros(numWinPerSub, numSeed);
    for i=1:numWinPerSub
        indxAllSeeds(i, :)=[indxSeed1(i), indxSeed2(i), indxSeed3(i), indxSeed4(i)];
        numUniqStat4Seeds(i,j)=length(unique(indxAllSeeds(i,:)));
    end
    for k=1:possibNumStatesAtATime
        prctDifNumStates(j,k)=length(find(numUniqStat4Seeds(:,j)==k))/numWinPerSub;
    end
end

close all
color=[255,69,0];
color=color/255;
figure(1)
num2plot=squeeze(numUniqStat4Seeds(:,j));
h=plot(Time,num2plot, '-','Color',color)
set(h,'LineWidth',lineWidth2)
ylim([0 5])
xlim([0 N_vol])
set(gca,'xTick', 0:200:800,'yTick', 0:5, 'XTickLabel',[],'YTickLabel',[],'box','off')
set(gca,'LineWidth',lineWidth1)
saveas(figure(1),[figDir, 'stateTC_sub', num2str(j), '_',session, 'sharedState.png'])