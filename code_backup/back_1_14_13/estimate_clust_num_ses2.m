clear
clc
close all

loadfile=load('/home/data/Projects/microstate/analysis_matlab/results/645/session2/z_Dyn_FC.mat')
        z_Dyn_FC=loadfile.z_Dyn_FC;
N_ROI=156;       
win_all_sub1=reshape(z_Dyn_FC(1,:,:),N_ROI,[])';
        win_all_sub2=reshape(z_Dyn_FC(2,:,:),N_ROI,[])';
        win_all_sub3=reshape(z_Dyn_FC(3,:,:),N_ROI,[])';
        win_all_sub4=reshape(z_Dyn_FC(4,:,:),N_ROI,[])';
        win_all_sub_seed=vertcat(win_all_sub1, win_all_sub2, win_all_sub3, win_all_sub4);
        disp ('Win of all seeds concatenated.')
save(['/home/data/Projects/microstate/analysis_matlab/results/645/session2/win_all_sub_seed.mat'],'win_all_sub_seed')

for N_clust=2:20
    disp(['Working on N_clust =', num2str(N_clust)])
        clear indx_stat ctrs_stat sumd D
        opts=statset('MaxIter',500);
        [indx_stat,ctrs_stat,sumd,D]=kmeans(win_all_sub_seed,N_clust,'Distance','correlation','emptyaction','singleton','Replicates',100,'options',opts);
        % compute the cluster validity index
         CVI(N_clust)=sum(sumd)./(sum(sum(D))-sum(sumd))
        index_stat(:,N_clust)=indx_stat;
                 %final_ctrs(N_ctrs,:)=ctrs_stat;
        %         final_ctrs_transp(:,:,n)=ctrs_stat';
        save(['/home/data/Projects/microstate/analysis_matlab/results/645/session2/indx_',N_clust,' clusters.mat'],'indx_stat')
        save(['/home/data/Projects/microstate/analysis_matlab/results/645/session2/state_',N_clust,' clusters.mat'],'ctrs_stat')
        disp (['k-means clustering with ',num2str(N_clust),' cluster done.'])
end
         figure(2)
         plot(CVI)
            grid on
            ylabel('CVI(Within/betw Dist)')
            xlabel('Number of Clusters')
            xlim([2 20])
            set(gca,'XTick',2:2:20,'YTick',0:0.2:1)
            title('Cluster Validity Index (CVI)')
            
       saveas(figure(1),['/home/data/Projects/microstate/analysis_matlab/fig/session2/CVI_all_win_seed.png']);