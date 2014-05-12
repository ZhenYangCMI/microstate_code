clear
clc

sessionList={'session1','session2'};
numSession=length(sessionList)
TR={'645'};

for i=1: numSession
    session=char(sessionList{i})
    tmp1=textread(['/home/data/Projects/microstate/subListLowMotion', session,'.txt'],'%s')
    tmp2=textread(['/home/data/Projects/microstate/subListHighMotion', session,'.txt'],'%s')
    subList=vertcat(tmp1, tmp2);
    numSub=length(subList)
    tmp3=textread(['/home/data/Projects/microstate/lambdaAtMaxLogLowMotion', session,'.txt'],'%f')
    tmp4=textread(['/home/data/Projects/microstate/lambdaAtMaxLogHighMotion', session,'.txt'],'%f')
    optimaLambdaList=vertcat(tmp3, tmp4);
    numoptimalLambda=length(optimaLambdaList)
    for j=1:numSub
        sub=char(subList(j))
        optimaLambda=optimaLambdaList(j)
        % the following function estimate the sensitivity of lambda to
        % paritial correlation windows using three similarity/distance
        % metrics. Type 3; compare the parital correlation matrix at
        % differen lambda with the partial correlation matrix computed at the lambdaAtMaxLog
        
        [correlEachSeedType3, distEachSeedType3,mutualInfoEachSeedType3]=lambdaSensitivAnalySeedSep(TR,session, sub, optimaLambda);
        
    end
end


% the code below is used to delete the unused files
analyDir=['/home/data/Projects/microstate/DPARSF_preprocessed/'];


for i=1: numSession
    session=char(sessionList{i})
    tmp1=textread(['/home/data/Projects/microstate/subListLowMotion', session,'.txt'],'%s')
    tmp2=textread(['/home/data/Projects/microstate/subListHighMotion', session,'.txt'],'%s')
    subList=vertcat(tmp1, tmp2);
    numSub=length(subList)
    
    for j=1:numSub
        sub=char(subList(j))
        similarityMetric={'correl','dist','mutualInfo'};
        numMetric=length(similarityMetric);
        compareType={'Type1','Type2','Type3'};
        numType=length(compareType);
        resultDir=[analyDir,'results/',char(TR), '/',session,'/lambdaSensitivity/seedFCWin/', sub,'/'];
        for m=1:numMetric
            metric=char(similarityMetric{m});
            for n=1:numType
                Type=char(compareType{n});
                delete([resultDir,metric,'AllSeeds',Type,'.mat'])
            end
        end             
        
    end
end


