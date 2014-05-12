function [seedPrct, subPrct, subSeedPrct] = DoQuickStats(indx)

numSubs = 22;
numStates = 7;
numWin = 128;
numSeeds = 4;

seedInd = [];
subInd = [];

for i = 1:numSeeds
    for j = 1: numSubs
        subInd = [subInd; ones([numWin 1]) * j];
    end
    
    seedInd = [seedInd; ones([numWin * numSubs 1]) * i];
end

for i = 1:numStates
    for j = 1:numSeeds
        seedPrct(i,j) = sum(indx(find(seedInd == j)) == i) / (numWin*numSubs); 
    end
    
    for j = 1:numSubs
        subPrct(i,j) = sum(indx(find(subInd == j)) == i) / (numWin*numSeeds); 
    end

end

subSeedPrct=zeros(numStates, numSeeds, numSubs);
for i = 1:numStates
    for j = 1:numSeeds
        for k = 1:numSubs
            begin_ndx=(j-1)*numSubs*numWin+(k-1)*numWin+1;
            end_ndx=begin_ndx+numWin-1;
            subSeedPrct(i,j,k)=sum(indx(begin_ndx:end_ndx)==i)/numWin;
            %disp(sprintf('begin %d, end %d, val %f',begin_ndx,end_ndx,subSeedPrct(i,j,k)));
        end
    end
end

        
            



