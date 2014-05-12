% this script compute the % time, all four seeds are in state 1-5
clear
clc
session='session1';
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/GSR/', session,'/'];
indx=load([resultDir, 'clustIndxNormWinAllSeeds_FullCorLasso_',session,'_10min.txt']);
numState=length(unique(indx));
figDir='/home2/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/GSR/session1/motion/';
numWinPerSeed=5984;

seed1=indx(1:numWinPerSeed);
seed2=indx(1+numWinPerSeed: numWinPerSeed*2);
seed3=indx(1+numWinPerSeed*2: numWinPerSeed*3);
seed4=indx(1+numWinPerSeed*3: numWinPerSeed*4);

numSeedsShareState=4;

numWinSeedsShareState=zeros(numState, 1);

for i=1:numState
    e=zeros(numWinPerSeed, 1);
    a=find(seed1==i);
    b=find(seed2==i);
    c=find(seed3==i);
    d=find(seed4==i);
    
   for j=1:length(a)
       e(a(j))=e(a(j))+1;
   end
   
   for k=1:length(b)
       e(b(k))=e(b(k))+1;
   end
   
   for m=1:length(c)
       e(c(m))=e(c(m))+1;
   end
   
   for n=1:length(d)
       e(d(n))=e(d(n))+1;
   end
   e(find(e==0))=[];
   
   %numIndx=length(find(e==numSeedsShareState)) %e.g percent time 3 seeds share a state
   numIndx=length(find(e>=numSeedsShareState)) %e.g percent time 3 or more seeds share a state
      numWinSeedsShareState(i, 1)=numIndx
end

prctInOneState=numWinSeedsShareState/numWinPerSeed*100

figure(1)
bar(prctInOneState)
xlabel('States')
title(['Percent time ', num2str(numSeedsShareState), 'seeds share the same state.'])
ylabel('%')
ylim([0 30])
%saveas(figure(1), [figDir, 'PercentTime_', num2str(numSeedsShareState),  'OrMoreSeedsShareState.png'])


lineWidth=2
    figure(2)
    h=bar(prctInOneState, 0.45) % TCPrctOverlapStateSum used 0.6, other 2 measures used 0.4 to make them the same width
    set(gca,'LineWidth',lineWidth)
        set(h,'LineWidth',lineWidth,'Facecolor',[1 0 0],'EdgeColor',[1 0 0])
   set(gca,'xTick',1:5, 'yTick', 0:5:20, 'box','off', 'xTickLabel',[],'yTickLabel',[]);
        ylim([0 20])
        saveas(figure(2), sprintf('/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/all_10min/GSR/%s/percent4seedsInOneState_%s.png', session, session))




