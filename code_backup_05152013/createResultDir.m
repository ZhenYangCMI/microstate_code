function createResultDir(RTList, sesseionList)

TRList={'645','2500'};
sessionList={'session1','session2'};

resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results'];
figDir= ['/home/data/Projects/microstate/DPARSF_preprocessed/fig'];
%analyDir=['/Users/zhenyang/Documents/microstate/'];

numTR=length(TRList);
numSession=length(sessionList);


for i=1:numTR
    TR=TRList{i};
    if ~exist([resultDir,'/',TR], 'dir')
        mkdir(resultDir,TR)
    end
    if ~exist([figDir,'/',TR], 'dir')
        mkdir(figDir,TR)
    end
    for j=1:numSession
        session=sessionList{j};
        if ~exist([resultDir,'/',TR,'/',session], 'dir')
            mkdir([resultDir,'/',TR], session)
        end
        if ~exist([figDir,'/',TR,'/',session], 'dir')
            mkdir([figDir,'/',TR], session)
        end
    end
end

