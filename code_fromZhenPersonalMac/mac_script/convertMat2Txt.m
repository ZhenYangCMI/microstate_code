clear
clc

session='session2';
numSub=22;
measureList1={'prctDifNumStates','transitions','TCPrctOverlapStateSum'};
measureList2={'subSeedPrct','stateFreq','duration','TCPrctOverlapEachState'};

for i=1:length(measureList1)
    measure=char(measureList1{i})
    inFileName=sprintf('%s_%s.mat',measure, session)
    tmp=load(inFileName);
    measure2D=tmp.(measure)
    outFileName=sprintf('%s_%s.txt',measure, session)
    save(outFileName,'-ascii', '-tabs','measure2D')
end

for i=1:length(measureList2)
    measure=char(measureList2{i})
    inFileName=sprintf('%s_%s.mat',measure, session)
    tmp=load(inFileName);
    measure3D=tmp.(measure)
    measure2D=reshape(measure3D,[],numSub)'
    outFileName=sprintf('%s_%s.txt',measure, session)
    save(outFileName,'-ascii', '-tabs','measure2D')
end