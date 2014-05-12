
clear
clc
close all

session='session1'
dataLength='all_10min';
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
% seedPairs=seed1 and seed2; 1&3; 1&4; 2&3; 2&4; 3&4 
numSeedPair=6;

if strcmp(session, 'session1')
    numState=5;
else
    numState=6;
end

numROI=156;
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,filesep,session,'/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep,session, filesep, 'stateTC/individual/'];

indx=load([resultDir,'clustIndxNormWinAllSeeds_FullCorLasso_',session,'_10min.txt']);
numClust=length(unique(indx));
disp ('Files loaded successfully.')

indxRecode=zeros(length(indx),1);

% match 5 pairs of states between two sessions and reorder them according
% to the strength of the correlation coefficient (i.e. new state 1-5 in
% session1 matches new state 1-5 in session2, and the correlation between
% them decrease from state 1 to state 5. State 6 in session2 is a session
% unique state)
if strcmp(session,'session1')
    indxRecode(find(indx==1))=4;
    indxRecode(find(indx==2))=2;
    indxRecode(find(indx==3))=5;
    indxRecode(find(indx==4))=1;
    indxRecode(find(indx==5))=3;
else
    indxRecode(find(indx==1))=1;
    indxRecode(find(indx==2))=2;
    indxRecode(find(indx==3))=4;
    indxRecode(find(indx==4))=6;
    indxRecode(find(indx==5))=3;
    indxRecode(find(indx==6))=5;
end

% compute for each subject and each seed: 1.the totoal number of
% transitions; 2. the frquency of each state
transitions=zeros(numSub,numSeed);
stateFreq=zeros(numState, numSeed, numSub);

for i=1:numSeed
    indxRecodeSeed=indxRecode((1+numWinPerSeed*(i-1)):numWinPerSeed*i);
    for j=1:numSub
        indxRecodeSeedSub=indxRecodeSeed((1+numWinPerSub*(j-1)):numWinPerSub*j);
        m=0;
        if strcmp(session,'session1')
            a(1:5)=0;
        else
            a(1:6)=0;
        end
        a(indxRecodeSeedSub(1))=a(indxRecodeSeedSub(1))+1;
        for k=1:numWinPerSub
            if k+1<=numWinPerSub
                if indxRecodeSeedSub(k+1)~=indxRecodeSeedSub(k)
                    m=m+1;
                    a(indxRecodeSeedSub(k+1))=a(indxRecodeSeedSub(k+1))+1;
                end
            else
                disp('k+1 exceeds the numWinPerSub.')
            end
        end
        transitions(j,i)=m;
        for n=1:numState
            stateFreq(n,i,j)=a(n);
        end
    end
end
save([resultDir, 'transitions_', session, '.mat'], 'transitions')
save([resultDir, 'stateFreq_', session, '.mat'], 'stateFreq')


% compute the %time spent in each state for each seed and each sub
subSeedPrct=zeros(numState, numSeed, numSub);
for n = 1:numState
    for i = 1:numSeed
        for j = 1:numSub
            begin_ndx=(i-1)*numSub*numWinPerSub+(j-1)*numWinPerSub+1;
            end_ndx=begin_ndx+numWinPerSub-1;
            subSeed(n,i,j)=sum(indxRecode(begin_ndx:end_ndx)==n);
            subSeedPrct(n,i,j)=sum(indxRecode(begin_ndx:end_ndx)==n)/numWinPerSub;
            %disp(sprintf('begin %d, end %d, val %f',begin_ndx,end_ndx,subSeedPrct(i,j,k)));
        end
    end
end
save([resultDir, 'subSeedPrct_', session, '.mat'], 'subSeedPrct')


% compute the microstate duration in unit of windows
for n = 1:numState
    for i = 1:numSeed
        for j = 1:numSub
            if stateFreq(n,i,j)~=0;
                duration(n,i,j)=subSeedPrct(n,i,j)/stateFreq(n,i,j)
            else
                duration(n,i,j)=0;
            end
        end
    end
end
save([resultDir, 'duration_', session, '.mat'], 'duration')

% comput the % overlap between each pair of seeds
prctOverlap=zeros(numSub, numState*numSeedPair)
for i=1:numSeed
    indxRecodeSeed=indxRecode((1+numWinPerSeed*(i-1)):numWinPerSeed*i);
    for j=1:numSub
        indxRecodeSeedSub=indxRecodeSeed((1+numWinPerSub*(j-1)):numWinPerSub*j);
        m=0;
        if strcmp(session,'session1')
            a(1:5)=0;
        else
            a(1:6)=0;
        end










