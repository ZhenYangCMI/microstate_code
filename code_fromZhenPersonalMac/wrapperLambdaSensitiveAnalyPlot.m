clear
clc
close all

sessionList={'session1','session2'};
numSession=length(sessionList)
TR={'645'};

for i=1:numSession
    session=char(sessionList{i})
    tmp1=textread(['/home/data/Projects/microstate/subListLowMotion', session,'.txt'],'%s')
    tmp2=textread(['/home/data/Projects/microstate/subListHighMotion', session,'.txt'],'%s')
    subList=vertcat(tmp1, tmp2);
    numSub=length(subList)
    for j=1:numSub
        sub=char(subList(j))
        disp(['Working on sub ', sub])
        [totoalNumFigs]=lambdaSensitivAnalyPlot(TR,session, sub)
    end
end

