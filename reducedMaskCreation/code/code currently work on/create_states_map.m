
indx=load('/home/data/Projects/microstate/analysis_matlab/results/645/session1/ses1_cor_all_seed/indx_all_sub_all_seeds.txt');
state=load('/home/data/Projects/microstate/analysis_matlab/results/645/session1/ses1_cor_all_seed/states_all_sub_all_seeds.txt')
temp=load('/home/data/Projects/microstate/analysis_matlab/results/645/session1/win_all_sub_seed.mat');
winSeedAll=temp.win_all_sub_seed;
numWinPerSub=128;
numSub=22;
numROI=156;
numState=7;
meanStatePerSub=zeros(numState,numROI,numSub);
for i=1:numSub
    winSub=winSeedAll((128*(i-1)+1):128*i,:);
    indxSub=indx(((128*(i-1)+1):128*i),1);
    for j=1:numState
        stateIndx=find(indxSub==j);
        stateMean= mean(winSub(stateIndx,:));
        for k=1:numROI
            meanStatePerSub(j:k:i)=stateMean(1,k);
        end
    end
end
save('/home/data/Projects/microstate/analysis_matlab/results/645/session1/ses1_cor_all_seed/avg_state_each_sub.mat','avg_state')

for j=1:numState
    stateAllSub=reshape(meanStatePerSub(j,:,:),numROI,[])';
    for k=1:numROI
        [~, pCurr,~, tstat]=ttest(stateAllSub(:,k)');
        t(j,k) = tstat.tstat;
        pCurr
    end
    
    %        pID(j,:) =
    %        temp = FDR(p(j,:),0.05);
end



