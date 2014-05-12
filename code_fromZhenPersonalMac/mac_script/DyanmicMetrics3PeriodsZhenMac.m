
clear
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

numWin1=ceil(numWinPerSub/3);
numWin2=numWinPerSub-2*ceil(numWinPerSub/3);
numWin3=ceil(numWinPerSub/3);

% comput the % overlap between each pair of seeds
seedPairList=[1, 2; 1, 3; 1, 4; 2, 3; 2,4; 3, 4];
numSeedPair=size(seedPairList,1);
TCPrctOverlap3Period1=zeros(numSub, numSeedPair);

disp('Period 1')
for i=1:numSeedPair
    seedPair=seedPairList(i,:);
    m=seedPair(1);
    n=seedPair(2);
    for j=1:numSub
        disp(['Working on seed pair', num2str(i), '_sub ', num2str(j)])
        begin_ndx1=(m-1)*numSub*numWinPerSub+(j-1)*numWinPerSub+1;
        end_ndx1=begin_ndx1+numWin1-1;
        begin_ndx2=(n-1)*numSub*numWinPerSub+(j-1)*numWinPerSub+1;
        end_ndx2=begin_ndx2+numWin1-1;
        disp(sprintf('begin1 %d, end1 %d',begin_ndx1,end_ndx1));
        disp(sprintf('begin2 %d, end2 %d',begin_ndx2,end_ndx2));
        indxFile1=indxRecode(begin_ndx1:end_ndx1);
        indxFile2=indxRecode(begin_ndx2:end_ndx2);
        t=length(find(indxFile1==indxFile2));
        TCPrctOverlap3Period1(j,i)=t/numWin1;
    end
end

disp('Period 2')
for i=1:numSeedPair
    seedPair=seedPairList(i,:);
    m=seedPair(1);
    n=seedPair(2);
    for j=1:numSub
        disp(['Working on seed pair', num2str(i), '_sub ', num2str(j)])
        begin_ndx1=(m-1)*numSub*numWinPerSub+(j-1)*numWinPerSub+numWin1+1;
        end_ndx1=begin_ndx1+numWin2-1;
        begin_ndx2=(n-1)*numSub*numWinPerSub+(j-1)*numWinPerSub+numWin1+1;
        end_ndx2=begin_ndx2+numWin2-1;
        disp(sprintf('begin1 %d, end1 %d',begin_ndx1,end_ndx1));
        disp(sprintf('begin2 %d, end2 %d',begin_ndx2,end_ndx2));
        indxFile1=indxRecode(begin_ndx1:end_ndx1);
        indxFile2=indxRecode(begin_ndx2:end_ndx2);
        t=length(find(indxFile1==indxFile2));
        TCPrctOverlap3Period2(j,i)=t/numWin2;
    end
end

disp('Period 3')
for i=1:numSeedPair
    seedPair=seedPairList(i,:);
    m=seedPair(1);
    n=seedPair(2);
    for j=1:numSub
        disp(['Working on seed pair', num2str(i), '_sub ', num2str(j)])
        begin_ndx1=(m-1)*numSub*numWinPerSub+(j-1)*numWinPerSub+numWin1+numWin2+1;
        end_ndx1=begin_ndx1+numWin3-1;
        begin_ndx2=(n-1)*numSub*numWinPerSub+(j-1)*numWinPerSub+numWin1+numWin2+1;
        end_ndx2=begin_ndx2+numWin3-1;
        disp(sprintf('begin1 %d, end1 %d',begin_ndx1,end_ndx1));
        disp(sprintf('begin2 %d, end2 %d',begin_ndx2,end_ndx2));
        indxFile1=indxRecode(begin_ndx1:end_ndx1);
        indxFile2=indxRecode(begin_ndx2:end_ndx2);
        t=length(find(indxFile1==indxFile2));
        TCPrctOverlap3Period3(j,i)=t/numWin3;
    end
end

save([resultDir, 'TCPrctOverlap3Period1_', session, '.mat'], 'TCPrctOverlap3Period1')
save([resultDir, 'TCPrctOverlap3Period2_', session, '.mat'], 'TCPrctOverlap3Period2')
save([resultDir, 'TCPrctOverlap3Period3_', session, '.mat'], 'TCPrctOverlap3Period3')




% comput the total number of transitions
transitions1=zeros(numSub,numSeed);
transitions2=zeros(numSub,numSeed);
transitions3=zeros(numSub,numSeed);

for i=1:numSeed
    for j=1:numSub
        disp('Period 1')
        disp(['Working on seed ', num2str(i), '_sub ', num2str(j)])
        begin_ndx=(i-1)*numSub*numWinPerSub+(j-1)*numWinPerSub+1;
        end_ndx=begin_ndx+numWin1-1;
        disp(sprintf('begin %d, end %d',begin_ndx,end_ndx));
        indxFile=indxRecode(begin_ndx:end_ndx);
        m=0;
        for k=1:numWin1
            if k+1<=numWin1
                if indxFile(k+1)~=indxFile(k)
                    m=m+1;
                end
            end
        end
        transitions1(j,i)=m;
    end
end

for i=1:numSeed
    for j=1:numSub
        disp('Period 2')
        disp(['Working on seed ', num2str(i), '_sub ', num2str(j)])
        begin_ndx=(i-1)*numSub*numWinPerSub+(j-1)*numWinPerSub+numWin1+1;
        end_ndx=begin_ndx+numWin2-1;
        disp(sprintf('begin %d, end %d',begin_ndx,end_ndx));
        indxFile=indxRecode(begin_ndx:end_ndx);
        m=0;
        for k=1:numWin2
            if k+1<=numWin2
                if indxFile(k+1)~=indxFile(k)
                    m=m+1;
                end
            end
        end
        transitions2(j,i)=m;
    end
end
for i=1:numSeed
    for j=1:numSub
        disp('Period 3')
        disp(['Working on seed ', num2str(i), '_sub ', num2str(j)])
        begin_ndx=(i-1)*numSub*numWinPerSub+(j-1)*numWinPerSub+numWin1+numWin2+1;
        end_ndx=begin_ndx+numWin3-1;
        disp(sprintf('begin %d, end %d',begin_ndx,end_ndx));
        indxFile=indxRecode(begin_ndx:end_ndx);
        m=0;
        for k=1:numWin3
            if k+1<=numWin3
                if indxFile(k+1)~=indxFile(k)
                    m=m+1;
                end
            end
        end
        transitions3(j,i)=m;
    end
end
transitions1
transitions2
transitions3
save([resultDir, 'transitions1_', session, '.mat'], 'transitions1')
save([resultDir, 'transitions2_', session, '.mat'], 'transitions2')
save([resultDir, 'transitions3_', session, '.mat'], 'transitions3')

% Divide the total time into 3 periods and compute the %time spent in each state for each seed and each sub for each period
subSeedPrct1=zeros(numState, numSeed, numSub);
subSeedPrct2=zeros(numState, numSeed, numSub);
subSeedPrct3=zeros(numState, numSeed, numSub);

% Period 1
for n = 1:numState
    for i = 1:numSeed
        for j = 1:numSub
            disp(['Working on state ', num2str(n), '_seed ', num2str(i), '_sub ', num2str(j)])
            begin_ndx=(i-1)*numSub*numWinPerSub+(j-1)*numWinPerSub+1;
            end_ndx=begin_ndx+numWin1-1;
            subSeedPrct1(n,i,j)=sum(indxRecode(begin_ndx:end_ndx)==n)/numWin1;
            disp(sprintf('begin %d, end %d, val %f',begin_ndx,end_ndx,subSeedPrct1(n,i,j)));
            
        end
    end
end

clear begin_ndx end_ndx
for n = 1:numState
    for i = 1:numSeed
        for j = 1:numSub
            disp(['Working on state ', num2str(n), '_seed ', num2str(i), '_sub ', num2str(j)])
            begin_ndx=(i-1)*numSub*numWinPerSub+(j-1)*numWinPerSub+numWin1+1;
            end_ndx=begin_ndx+numWin2-1;
            subSeedPrct2(n,i,j)=sum(indxRecode(begin_ndx:end_ndx)==n)/numWin2;
            disp(sprintf('begin %d, end %d, val %f',begin_ndx,end_ndx,subSeedPrct2(n,i,j)));
        end
    end
end

clear begin_ndx end_ndx
for n = 1:numState
    for i = 1:numSeed
        for j = 1:numSub
            disp(['Working on state ', num2str(n), '_seed ', num2str(i), '_sub ', num2str(j)])
            begin_ndx=(i-1)*numSub*numWinPerSub+(j-1)*numWinPerSub+numWin1+numWin2+1;
            end_ndx=begin_ndx+numWin3-1;
            subSeedPrct3(n,i,j)=sum(indxRecode(begin_ndx:end_ndx)==n)/numWin3;
            disp(sprintf('begin %d, end %d, val %f',begin_ndx,end_ndx,subSeedPrct3(n,i,j)));
        end
    end
end
save([resultDir, 'subSeedPrct1_', session, '.mat'], 'subSeedPrct1')
save([resultDir, 'subSeedPrct2_', session, '.mat'], 'subSeedPrct2')
save([resultDir, 'subSeedPrct3_', session, '.mat'], 'subSeedPrct3')


