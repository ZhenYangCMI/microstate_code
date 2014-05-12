

sessionList={'FunImgR','S2_FunImgR'};
numSession=length(sessionList);

clear
clc
dataDir=['/home/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/S2_FunImgR/'];
subList={'0021002', '0021006', '0021018', '0021024', '1427581', '1793622',...
    '1961098', '2475376', '2799329', '2842950', '3201815', '3313349'...
    '3315657', '3795193', '3808535', '3893245', '4176156', '4288245',...
    '7055197', '8574662', '8735778', '9630905'};
numSub=length(subList);

for i=1:numSub
    subDir=[dataDir, subList{i}];
    disp (subDir)
    delete ([subDir,'/rest.nii'],[subDir,'/rrest.mat'])
end


for i=1:numSession
    sessionDir=[dataDir, '/', sessionList{i}];
    for j=1:numSub
        subDir=[sessionDir,'/', subList{j}];
        disp (subDir)
        delete ([subDir,'/rrest.mat'],[subDir,'/rrest.nii'])
    end
end

