% this script compute the % time, all four seeds are in state 1-5
clear
clc
session='session1';
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/all_10min/GSR/', session,'/'];
indx=load([resultDir, 'clustIndxNormWinAllSeeds_FullCorLasso_',session,'_10min.txt']);
numState=length(unique(indx));

numWinPerSeed=5984;

seed1=indx(1:numWinPerSeed);
seed2=indx(1+numWinPerSeed: numWinPerSeed*2);
seed3=indx(1+numWinPerSeed*2: numWinPerSeed*3);
seed4=indx(1+numWinPerSeed*3: numWinPerSeed*4);

prctInOneState=zeros(numState, 1);
for i=1:numState
    a=find(seed1==i);
    b=find(seed2==i);
    c=find(seed3==i);
    d=find(seed4==i);
    t=0;
    for j=1:length(a)
        
        for k=1:length(b)
            if a(j)~=b(k)
                continue
            end
            for m=1:length(c)
                if a(j)~=c(m)
                    continue
                end
                for n=1:length(d)
                    if a(j)~=d(n)
                        continue
                    else
                        a(j)
                        b(k)
                        c(m)
                        d(n)
                        t=t+1
                    end
                end
            end
        end
    end
    prctInOneState(i, 1)=t
end

prctInOneState=prctInOneState/numWinPerSeed*100
bar(prctInOneState)
xlabel('States')
title('Percent time all four seeds in the same state.')
ylabel('%')




