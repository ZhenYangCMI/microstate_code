
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

if strcmp(session, 'session1')
    numState=5;
else
    numState=6;
end

numROI=156;
resultDir=['/Users/zhenyang/Desktop/Zhen/results_4_3_13/all_10min/',filesep,session,'/'];
figDir=['/Users/zhenyang/Desktop/Zhen/figs/'];

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

% compute the transition matrix of Markov Chain
TM=zeros(numState,numState,numSub,numSeed);

for i=1:numSeed
    indxRecodeSeed=indxRecode((1+numWinPerSeed*(i-1)):numWinPerSeed*i);
    for j=1:numSub
        indxRecodeSeedSub=indxRecodeSeed((1+numWinPerSub*(j-1)):numWinPerSub*j);
        markovChain=indxRecodeSeedSub;
        Norder=1;
        [ transitionMatrix columnStates ] = getTransitionMatrix(markovChain,Norder);
        if strcmp(session,'session1')
            maxIndx=5;
        else
            maxIndx=6;
        end
        if max(indxRecodeSeedSub)==maxIndx
            TM(:,:,j,i)=transitionMatrix'; %the row represent state at time t and each column the time t+1
        elseif max(indxRecodeSeedSub)==maxIndx-1
            x=zeros(1,maxIndx-1);
            y=zeros(maxIndx,1);
            transitionMatrix=vertcat(transitionMatrix,x);
            transitionMatrix=horzcat(transitionMatrix,y);
            TM(:,:,j,i)=transitionMatrix';
        elseif max(indxRecodeSeedSub)==maxIndx-2
            x=zeros(2,maxIndx-2);
            y=zeros(maxIndx,2);
            transitionMatrix=vertcat(transitionMatrix,x);
            transitionMatrix=horzcat(transitionMatrix,y);
            TM(:,:,j,i)=transitionMatrix';
        elseif max(indxRecodeSeedSub)==maxIndx-3
            x=zeros(3,maxIndx-3);
            y=zeros(maxIndx,3);
            transitionMatrix=vertcat(transitionMatrix,x);
            transitionMatrix=horzcat(transitionMatrix,y);
            TM(:,:,j,i)=transitionMatrix';
        elseif max(indxRecodeSeedSub)==maxIndx-4
            x=zeros(4,maxIndx-4);
            y=zeros(maxIndx,4);
            transitionMatrix=vertcat(transitionMatrix,x);
            transitionMatrix=horzcat(transitionMatrix,y);
            TM(:,:,j,i)=transitionMatrix';
            
        end
    end
end

% compute the probablity of the transition matrix
probTM=zeros(numState,numState,numSub,numSeed);
for i=1:numSeed
    for j=1:numSub
        rowSumTM= sum(TM(:,:,j,i)')';
        probTM(:,:,j,i)=TM(:,:,j,i)./repmat(rowSumTM,1,maxIndx);
        probTM(isnan(probTM))=0;
    end
end

meanTMSeed=zeros(numState,numState,numSeed);
for i=1:numSeed
    TMSeed=squeeze(probTM(:,:,:,i));
    TMSeed2D=reshape(TMSeed,[],numSub)';
    tmp=mean(TMSeed2D);
    meanProbTM(:,:,i)=reshape(tmp,numState,numState);
end

meanProbTMDifState=zeros(numState,numState,numSeed);
for i=1:numSeed
    data=squeeze(meanProbTM(:,:,i));
    % remove the probability of entering or exiting its own state
    for j=1:numState
        for k=1:numState
            if j==k
                data(j,k)=0;
            end
        end
    end
    for j=1:numState
        for k=1:numState
            if j~=k
                meanProbTMDifState(j,k,i)=data(j,k)/sum(data(j,:));
            end
        end
    end
end
meanProbTMDifState

save([resultDir, 'meanProbTM_', session, '.mat'], 'meanProbTM')
save([resultDir, 'probTM_', session, '.mat'], 'probTM')
save([resultDir, 'meanProbTMDifState_', session, '.mat'], 'meanProbTMDifState')

%        for i=1:numSeed
%         [Eig_Vector,Eig_value] = eig(meanTM(:,:,i));
%        end
%         stat_vector=Eig_Vector(:,1)';
%
%         figure(17+4*(i-1))
%         plot(stat_vector,'-or','Markersize', 10)
%         ylabel('Stationary Probability')
%         xlabel('States')
%         title('Steady-state Behavior')



% comput the % overlap between each pair of seeds
seedPairList=[1, 2; 1, 3; 1, 4; 2, 3; 2,4; 3, 4];
numSeedPair=size(seedPairList,1);
TCPrctOverlapEachState=zeros(numState, numSeedPair, numSub);

for i=1:numSeedPair
    seedPair=seedPairList(i,:);
    m=seedPair(1);
    n=seedPair(2);
    seed1=indxRecode((1+numWinPerSeed*(m-1)):numWinPerSeed*m);
    seed2=indxRecode((1+numWinPerSeed*(n-1)):numWinPerSeed*n);
    for j=1:numSub
        seedSub1=seed1((1+numWinPerSub*(j-1)):numWinPerSub*j);
        seedSub2=seed2((1+numWinPerSub*(j-1)):numWinPerSub*j);
        for k=1:numState
            t=0;
            tmp1=find(seedSub1==k);
            tmp2=find(seedSub2==k);
            if isempty(tmp1);
                t=0;
            else
                for p=1:length(tmp1)
                    if ismember(tmp1(p),tmp2(:))
                        tmp2Overlap=tmp2(find(tmp2==tmp1(p)));
                        disp(['seed', num2str(m), 'indx ', num2str(tmp1(p)) ' overlapped with seed', num2str(n), 'indx ', num2str(tmp2Overlap)])
                        t=t+1;
                    end
                end
            end
            TCPrctOverlapEachState(k,i,j)=t/numWinPerSub;
        end
    end
end
TCPrctOverlapStateSum=squeeze(sum(TCPrctOverlapEachState))';
save([resultDir, 'TCPrctOverlapEachState_', session, '.mat'], 'TCPrctOverlapEachState')
save([resultDir, 'TCPrctOverlapStateSum_', session, '.mat'], 'TCPrctOverlapStateSum')



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
            probTM(1:5)=0;
        else
            probTM(1:6)=0;
        end
        probTM(indxRecodeSeedSub(1))=probTM(indxRecodeSeedSub(1))+1;
        for k=1:numWinPerSub
            if k+1<=numWinPerSub
                if indxRecodeSeedSub(k+1)~=indxRecodeSeedSub(k)
                    m=m+1;
                    probTM(indxRecodeSeedSub(k+1))=probTM(indxRecodeSeedSub(k+1))+1;
                end
            else
                disp('k+1 exceeds the numWinPerSub.')
            end
        end
        transitions(j,i)=m;
        for n=1:numState
            stateFreq(n,i,j)=probTM(n);
        end
    end
end
% save([resultDir, 'transitions_', session, '.mat'], 'transitions')
% save([resultDir, 'stateFreq_', session, '.mat'], 'stateFreq')


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
% save([resultDir, 'subSeedPrct_', session, '.mat'], 'subSeedPrct')


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










