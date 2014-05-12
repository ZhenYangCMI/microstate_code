clear
clc
close all

loadfile=load('/home/data/Projects/microstate/analysis_matlab/results/645/session1/z_Dyn_FC.mat')
zDynFC=loadfile.z_Dyn_FC;
numROI=156;
winSeed1=reshape(zDynFC(1,:,:),numROI,[])';
winSeed2=reshape(zDynFC(2,:,:),numROI,[])';
winSeed3=reshape(zDynFC(3,:,:),numROI,[])';
winSeed4=reshape(zDynFC(4,:,:),numROI,[])';
winSeedAll=vertcat(winSeed1, winSeed2, winSeed3, winSeed4);
disp ('Win of all seeds concatenated.')
save(['/home/data/Projects/microstate/analysis_matlab/results/645/session1/win_all_sub_seed.mat'],'win_all_sub_seed')
for numClust=2:20
    disp(['Working on nuCclust =', num2str(numClust)])
    clear stateIndx stateCtrs sumd D
    opts=statset('MaxIter',500);
    [stateIndx,stateCtrs,sumd,D]=kmeans(winSeedAll,numClust,'Distance','correlation','emptyaction','singleton','Replicates',100,'options',opts);
    % compute the cluster validity index
    CVI(numClust)=sum(sumd)./(sum(sum(D))-sum(sumd))
    stateIndx(:,numClust)=stateIndx;
    %final_ctrs(:,:,N_clust)=ctrs_stat;
    %         final_ctrs_transp(:,:,n)=ctrs_stat';
    save(['/home/data/Projects/microstate/analysis_matlab/results/645/session1/indx_',num2str(numClust),' clusters.mat'],'stateIndx')
    save(['/home/data/Projects/microstate/analysis_matlab/results/645/session1/state_',num2str(numClust),' clusters.mat'],'stateCtrs')
    disp (['k-means clustering with',num2str(numClust),' cluster done.'])
end
figure(1)
plot(CVI)
grid on
ylabel('CVI(Within/betw Dist)')
xlabel('Number of Clusters')
xlim([2 20])
set(gca,'XTick',2:2:20,'YTick',0:0.2:1)
title('Cluster Validity Index (CVI)')

saveas(figure(1),['/home/data/Projects/microstate/analysis_matlab/fig/session1/CVI_all_win_seed.png']);