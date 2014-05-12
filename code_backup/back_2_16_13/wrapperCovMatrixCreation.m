clear
clc

sessionList={'session1','session2'};
numSession=length(sessionList)
TR={'645'};

for i=1: numSession
    session=char(sessionList{i})
    tmp1=textread(['/home/data/Projects/microstate/subListLowMotion', session,'.txt'],'%s')
    tmp2=textread(['/home/data/Projects/microstate/subListHighMotion', session,'.txt'],'%s')
    subList=[tmp1, tmp2];
    numSub=length(subList)
    for j=1:numSub
        sub=char(subList(j))
        % the following function create and save the full correlation and full and  patiral
        % covariance matrix
        [fullCorrel, fullLasso, partialLasso]=covMatrixCreationSingleSub(TR,session,sub);
    end
end
