clear
clc

TR={'645'};
sessionList={'session1','session2'};
numSession=length(sessionList);


analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];

for i=1:numSession
    session=sessionList{i};
    resultDir1=[analyDir,'results/',char(TR),'/',char(session),'/lambdaSensitivity/covMatrix'];
    resultDir2=[analyDir,'results/',char(TR), '/',char(session),'/lambdaSensitivity/seedFCWin'];
    tmp1=textread(['/home/data/Projects/microstate/subListLowMotion', char(session),'.txt'],'%s');
    tmp2=textread(['/home/data/Projects/microstate/subListHighMotion', char(session),'.txt'],'%s');
    subList=vertcat(tmp1, tmp2)
    numSub=length(subList)
    
    for j=1:numSub
        sub=char(subList(j));
        disp(['Work on ', session,' sub ', sub])
        mkdir(resultDir1,sub)
        mkdir(resultDir2,sub)
    end
end
