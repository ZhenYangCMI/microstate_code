clear
clc
close all

% window width in TRs
W_width = 69;
% window sliding step in TRs
step=3;
% time
T=0:W_width-1;
% repetition time in secs
TR=0.645;
% sampling rate
fs=1/TR;
% total number of time points
N_vol=450;
% number of windows
N_win=floor((N_vol-W_width)/step)+1;


% creat the rectangular window
RW = ones(1,W_width);
% figure(1); subplot(1,2,1); plot(T,RW)
% xlim([0 length(RW)-1])

% Calculate gaussian
% Gaussian Mean
Gm=(W_width-1)/2;
% Gaussian Standard Deviation in TRs
Gs=3;
GW=exp(-(T-Gm).^2/(2*Gs.^2));
% figure(2); subplot(1,2,1); plot(GW)
% xlim([0 length(GW)-1])

% convolve rectangular and gaussian
WC=conv(RW,GW);
% figure(3); subplot(1,2,1); plot(WC)
% xlim([0 length(WC)-1])

SubList={'0021002', '0021006', '0021018', '0021024', '1427581', '1793622',...
    '1961098', '2475376', '2799329', '2842950', '3201815', '3313349'...
    '3315657', '3795193', '3808535', '3893245', '4176156', '4288245',...
    '7055197', '8574662', '8735778', '9630905'};
SesList={'session1','session2'};
%TRList={'645','2500'};

% SubList={'0021002', '0021006', '0021018', '0021024', '1427581'};
% SesList={'session1'};
TRList={'645'};

N_sub=size(SubList,2);
N_ses=size(SesList,2);
N_TR=size(TRList,2);


ROI_index=load('/home/data/Projects/microstate/analysis_matlab/doc/Crad179_ROI_index.txt');
% rearrange the regions by lobes, after rearrange:lobe1=F (75ROIs, col 1:75); lobe2=T (31ROIs, col 76:106); lobe3=p(20ROIs,col 107:126); lobe4=o(21ROIs,col 127:147);
% lobe5=subcortical(9ROIs, col 148:156)
col_order=xlsread('/home/data/Projects/microstate/analysis_matlab/doc/Crad_179_lob_hier_label_final.xls','Crad179_ROI_index','a2:c157');
col_order=col_order(:,3)';
analysisdir=['/home/data/Projects/microstate/analysis_matlab/'];
figdir=['/home/data/Projects/microstate/analysis_matlab/fig/'];
maskdir=['/home/data/Projects/microstate/analysis_matlab/mask/']

for i=1:N_TR
    TR=TRList{i};
    TRdir=['/home/data/Projects/microstate/analysis_matlab/data/',TR];
    if ~exist(['/home/data/Projects/microstate/analysis_matlab/results/',TR], 'dir')
        mkdir('/home/data/Projects/microstate/analysis_matlab/results', TR)
    end
    TR_result_dir=['/home/data/Projects/microstate/analysis_matlab/results/',TR];
    for j=1:N_ses
        ses=SesList{j};
        sesdir=[TRdir,'/',ses];
        if ~exist([TR_result_dir,'/',ses], 'dir')
            mkdir(TR_result_dir, ses)
        end
        ses_result_dir=[TR_result_dir,'/',ses];
        
        for k=1:N_sub
            subdir=[sesdir,'/', char(SubList{k})];
            disp (['Working on sub ', char(SubList{k}),' ......'])
            
            seed_ROISignals = load([subdir,'/ROISignals_seed_ROISignal.mat']);
            TC1=seed_ROISignals.ROISignals;
            N_seed=size(TC1,2);
            
            ROI_ROISignals=load([subdir,'/ROISignals_atlas_ROISignal.mat']);
            TC2=ROI_ROISignals.ROISignals;
            N_ROI=size(TC2,2);
            
            % concatenate the time series of seeds and ROIs
            TC=[TC1,TC2];
            
            % Stationary FC estimation
            
            % Compute the full correlation matrix
            Corr_stat=zeros(size(TC,2),size(TC,2),N_sub);
            Corr_stat(:,:,k)=corrcoef(TC);
            
            % Extract the stationary FC between all pairs of seeds
            for n=1:N_seed
                for m=1:N_seed
                    stat_FC_seed(n,m,k)=Corr_stat(n,m,k);
                    z_stat_FC_seed(n,m,k)=0.5*log((1+stat_FC_seed(n,m,k))./(1-stat_FC_seed(n,m,k)));
                end
            end
            
            % Extract the stationary FC between seeds and ROIs
            stat_FC=zeros(N_seed,N_ROI,N_sub);
            for n=1:N_seed
                for m=N_seed+1:N_ROI+N_seed
                    stat_FC(n,m-N_seed,k)=Corr_stat(n,m,k);
                    z_stat_FC(n,m-N_seed,k)=0.5*log((1+stat_FC(n,m-N_seed,k))./(1-stat_FC(n,m-N_seed,k)));
                end
            end
            
        end
        
        % average across all subjects
        temp1=reshape(z_stat_FC_seed,[],k)';
        avg_temp1=mean(temp1);
        avg_z_stat_FC_seed=reshape(avg_temp1,N_seed,N_seed);
        save([ses_result_dir,'/avg_z_stat_FC_seed.mat'],'avg_z_stat_FC_seed');
        
        temp2=reshape(z_stat_FC,[],k)';
        avg_temp2=mean(temp2);
        avg_z_stat_FC=reshape(avg_temp2,N_seed,N_ROI);
        save([ses_result_dir,'/avg_z_stat_FC.mat'],'avg_z_stat_FC');
        
        % plot the stat_FC between seeds
        figure(1)
        imagesc(avg_z_stat_FC_seed)
        colorbar('EastOutside')
        set(gca,'XTick',1:4,'YTick',1:4);
        ylabel('Seeds')
        xlabel('Seeds')
        %xlim([-1 1]);
        title('Stat FC between Seeds (Z)')
        saveas(figure(1),[figdir, 'stat_FC_between_seeds.png']);
        
        % rearrange stat_FC ROIs by lobe
        
        avg_z_stat_FC_bylobe=avg_z_stat_FC(:,[col_order]);
        
        % Plot the stat_FC(Z) between seeds and ROIs
        figure(2)
        imagesc(avg_z_stat_FC_bylobe)
        colorbar('EastOutside')
        caxis([-0.5,0.5])
        set(gca,'YTickLabel','');
        ylabel('Seeds')
        xlabel('ROIs')
        title('Stat FC between Seeds and ROIs (Z)')
        saveas(figure(2),[figdir, 'stat_FC_seeds_and_ROIs.png']);
        
        disp('Stat_FC estimation done!')
        
        %% Dynamic FC estimation
        
        for k=1:N_sub
            subdir=[sesdir,'/', char(SubList{k})];
            disp (['Working on sub ', char(SubList{k}),' ......'])
            
            seed_ROISignals = load([subdir,'/ROISignals_seed_ROISignal.mat']);
            TC1=seed_ROISignals.ROISignals;
            N_seed=size(TC1,2);
            
            ROI_ROISignals=load([subdir,'/ROISignals_atlas_ROISignal.mat']);
            TC2=ROI_ROISignals.ROISignals;
            N_ROI=size(TC2,2);
            
            % concatenate the time series of seeds and ROIs
            TC=[TC1,TC2];
            
            
            % apply the sliding window to the time series
            asize = size(TC);
            WC_final=WC(1,ceil(W_width/2):(W_width+ceil(W_width/2)-1))';
            % win(W_width,asize(2),N_win_tot);
            for q=1+N_win*(k-1):N_win*k;
                for n=((q-1)-(k-1)*N_win)*step+1:W_width+((q-1)-(k-1)*N_win)*step
                    for m=1:asize(2)
                        win(n-((q-1)-(k-1)*N_win)*step,m,q)=TC(n,m)*WC_final(n-((q-1)-(k-1)*N_win)*step);
                    end
                end
            end
        end
        save([ses_result_dir,'/win_all_sub.mat'],'win');
        
        disp ('Window applying done!')
        
        % generate the full correlation matrix for each win and extrac the dyn FC
        
        % Dyn_FC(N_seed,N_ROI,N_win_tot);
        % Dyn_FC_seed(N_seed,N_seed,N_win_tot);
        N_win_tot=N_win*N_sub;
        for q=1:N_win_tot
            Corr_dyn(:,:,q)=corrcoef(win(:,:,q));
            for n=1:N_seed
                for m=1:N_seed
                    Dyn_FC_seed(n,m,q)=Corr_dyn(n,m,q);
                    z_Dyn_FC_seed(n,m,q)=0.5*log((1+Dyn_FC_seed(n,m,q))./(1-Dyn_FC_seed(n,m,q)));
                end
                for m=N_seed+1:N_ROI+N_seed
                    Dyn_FC(n,m-N_seed,q)=Corr_dyn(n,m,q);
                    z_Dyn_FC(n,m-N_seed,q)=0.5*log((1+Dyn_FC(n,m-N_seed,q))./(1-Dyn_FC(n,m-N_seed,q)));
                end
            end
        end
        save([ses_result_dir,'/Dyn_FC.mat'],'Dyn_FC');
        save([ses_result_dir,'/z_Dyn_FC.mat'],'z_Dyn_FC');
        disp ('Full correlation of each window computed!')
    end
end


% find the optimal number of clusters for each seed
for k=1:N_sub
    for n=1:N_seed
        win_sig_sub=reshape(z_Dyn_FC(n,:,(N_win*(k-1)+1):k*N_win),N_ROI,[])';
        clear indx ctrs sumd D
        % k-mean clustering
        for N_clust=2:10
            [indx,ctrs,sumd,D]=kmeans(win_sig_sub,N_clust,'Distance','cityblock','emptyaction','singleton','Replicates',500);
            
            % compute the cluster validity index
            CVI(n,N_clust,k)=sum(sumd)./(sum(sum(D))-sum(sumd));
        end
        opt_N_clust(k,n)=find(CVI(n,2:10,k)<0.05, 1);
        
        figure(8+k)
        subplot(2,2,n)
        plot(CVI(n,:,k))
        grid on
        ylabel('CVI(Within/betw Dist)')
        xlabel('Number of Clusters')
        xlim([2 10])
        set(gca,'XTick',2:10,'YTick',0:0.2:1)
        title('Cluster Validity Index (CVI)')
        
    end
    saveas(figure(8+k),[figdir, 'CVI_sub', num2str(k),'.png']);
end
save([ses_result_dir,'/opt_N_clust.txt'],'opt_N_clust','-ascii', '-double', '-tabs');

% Common number of state was defined as the most frequenty
% appreared solutions across all subjects and all seeds

N_stat=mode(mode(opt_N_clust));


for n=1:N_seed
    win_all_sub=reshape(z_Dyn_FC(n,:,:),N_ROI,[])';
    for N_clust=2:10
    [indx,ctrs,sumd,D]=kmeans(win_all_sub,N_clust,'Distance','cityblock','emptyaction','singleton','Replicates',500);
     % compute the cluster validity index
            CVI(n,N_clust,k)=sum(sumd)./(sum(sum(D))-sum(sumd));
    index_stat(:,n)=indx_stat;
    final_ctrs(:,:,n)=ctrs_stat;
    final_ctrs_transp(:,:,n)=ctrs_stat';
    
    final_ctrs_cat=reshape(final_ctrs_transp, N_ROI, []);
    
    % compare the similarity between states
    [stat_sim,p_val]=corrcoef(final_ctrs_cat);
end
save([analysisdir,'results/state_sim.txt'],'stat_sim','-ascii', '-double', '-tabs')
save([analysisdir,'results/p_val.txt'],'p_val','-ascii', '-double', '-tabs')


% plot the correlation of states between seeds
figure(31)
subplot(3,2,1)
imagesc(stat_sim(N_stat+1:N_stat*2,1:N_stat))
colorbar('EastOutside')
caxis([-1,1])
set(gca,'XTick',1:N_stat,'YTick', 1:N_stat);
ylabel('Seed 2 States')
xlabel('Seed 1 States')
title('State Similarity')

subplot(3,2,2)
imagesc(stat_sim(N_stat*2+1:N_stat*3,1:N_stat))
colorbar('EastOutside')
caxis([-1,1])
set(gca,'XTick',1:N_stat,'YTick', 1:N_stat);
ylabel('Seed 3 States')
xlabel('Seed 1 States')
title('State Similarity')

subplot(3,2,3)
imagesc(stat_sim(N_stat*3+1:N_stat*4,1:N_stat))
colorbar('EastOutside')
caxis([-1,1])
set(gca,'XTick',1:N_stat,'YTick', 1:N_stat);
ylabel('Seed 4 States')
xlabel('Seed 1 States')
title('State Similarity')

subplot(3,2,4)
imagesc(stat_sim(N_stat*2+1:N_stat*3,N_stat+1:N_stat*2))
colorbar('EastOutside')
caxis([-1,1])
set(gca,'XTick',1:N_stat,'YTick', 1:N_stat);
ylabel('Seed 3 States')
xlabel('Seed 2 States')
title('State Similarity')

subplot(3,2,5)
imagesc(stat_sim(N_stat*3+1:N_stat*4,N_stat+1:N_stat*2))
colorbar('EastOutside')
caxis([-1,1])
set(gca,'XTick',1:N_stat,'YTick', 1:N_stat);
ylabel('Seed 4 States')
xlabel('Seed 2 States')
title('State Similarity')

subplot(3,2,6)
imagesc(stat_sim(N_stat*3+1:N_stat*4,N_stat*2+1:N_stat*3))
colorbar('EastOutside')
caxis([-1,1])
set(gca,'XTick',1:N_stat,'YTick', 1:N_stat);
ylabel('Seed 4 States')
xlabel('Seed 3 States')
title('State Similarity')
saveas(figure(31),[figdir,'corr_stat_between_seeds.png']);

% plot the state
figure (32)
for n=1:N_seed
    subplot(2,2,n)
    final_ctrs_bylobe=final_ctrs(:,[col_order],:);
    imagesc(final_ctrs_bylobe(:,:,n))
    colorbar('EastOutside')
    ylabel('State')
    xlabel('ROIs')
    title(['Temporal Dynamics of FC for Seed',num2str(n)])
end
saveas(figure(32),[figdir,'states.png']);

[Outdata,VoxDim,Header]=rest_readfile([maskdir,'final_reduced.nii']);
[nDim1 nDim2 nDim3]=size(Outdata);
temp=unique(Outdata);
ROI_index=temp(find(temp~=0));
for n=1:N_seed
    for t=1:N_stat
        Statemap=Outdata;
        state_ctrs=final_ctrs(t,:,n);
        for m=1:N_ROI
            Statemap(find(Outdata==ROI_index(m)))=state_ctrs(m);
        end
        
        Header.pinfo = [1;0;0];
        Header.dt    =[16,0];
        rest_WriteNiftiImage(Statemap,Header,[figdir,'seed',num2str(n),'_state.nii_', num2str(t)]);
    end
end
end
end
% plot dyn FC exemplars, TC, PS for one sub
sub_indx=2;

% Plot examplar dynamic FC between seeds and ROIs in Z
z_Dyn_FC_bylobe=z_Dyn_FC(:,[col_order],:);
exemp_z_Dyn_FC=z_Dyn_FC(:,:,((sub_indx-1)*N_win+1):sub_indx*N_win);
exemp_z_Dyn_FC_bylobe=z_Dyn_FC_bylobe(:,:,((sub_indx-1)*N_win+1):sub_indx*N_win);
for q=1:20:N_win
    figure(3)
    subplot(3,3,(q-1)/20+1)
    imagesc(exemp_z_Dyn_FC_bylobe(:,:,q))
    colorbar('EastOutside')
    set(gca,'YTick',1:4);
    ylabel('Seeds')
    xlabel('ROIs')
    title('Dynamic FC (Z)')
end
saveas(figure(3),[figdir,'dyn_FC_seeds_ROIs.png']);

% Examplar time course and power sepctral of dynamic FC between seeds and ROIs

plotcolor=['r' 'g' 'b' 'c' 'm' 'k' 'y'];

for n=1:N_seed
    for m=1:50:N_ROI;
        FC_TC=reshape(exemp_z_Dyn_FC(n,m,:),1,[]);
        PS= abs(fftshift(fft(FC_TC)));
        figure(4)
        subplot(2,2,n)
        Time=ceil(W_width/2):step:(ceil(W_width/2)+step*(N_win-1));
        plot(Time,FC_TC,plotcolor(rem(m,size(plotcolor,2))+1))
        xlim([0,N_vol])
        ylabel('Correltion (Z)')
        xlabel('Time (TR)')
        title('Time Course of Dynamic FC')
        hold on
        figure(5)
        subplot(2,2,n)
        t=-fs/2:fs/(length(PS)):fs/2-fs/(length(PS));
        plot(t, PS, plotcolor(rem(m,size(plotcolor,2))+1))
        ylabel('Amplitude')
        xlabel('Frequency(Hz)')
        title('Power Spectral of Dynamic FC')
        hold on
    end
end
saveas(figure(4),[figdir,'TC_dyn_FC_seeds_ROIs.png']);
saveas(figure(5),[figdir,'PS_dyn_FC_seeds_ROIs.png']);

% Plot examplar dynamic FC between seeds in Z
exemp_z_Dyn_FC_seed=z_Dyn_FC_seed(:,:,((sub_indx-1)*N_win+1):sub_indx*N_win);
for q=1:20:N_win
    figure(6)
    subplot(3,3,(q-1)/20+1)
    imagesc(exemp_z_Dyn_FC_seed(:,:,q))
    colorbar('EastOutside')
    set(gca,'XTick',1:4,'YTick',1:4);
    ylabel('Seeds')
    xlabel('Seeds')
    title('Dynamic FC (Z)')
end
saveas(figure(6),[figdir,'dyn_FC_between_seeds.png']);

% Examplar time course and power sepctral of dynamic FC between seeds

plotcolor=['r' 'b' 'y' 'c' 'm' 'g' 'k'];

for n=1:N_seed-1
    for m=n+1:N_seed;
        FC_seed_TC=reshape(exemp_z_Dyn_FC_seed(n,m,:),1,[]);
        PS= abs(fftshift(fft(FC_seed_TC)));
        figure(7)
        Time=ceil(W_width/2):step:(ceil(W_width/2)+step*(N_win-1));
        plot(Time,FC_seed_TC,plotcolor(rem(2*n+m,size(plotcolor,2))+1))
        xlim([0,N_vol])
        ylabel('Correltion (Z)')
        xlabel('Time (TR)')
        title('Time Course of Dynamic FC between Seeds')
        legend('S1-S2','S1-S3','S1-S4','S2-S3','S2-S4','S3-S4')
        hold on
        
        figure(8)
        t=-fs/2:fs/(length(PS)):fs/2-fs/(length(PS));
        plot(t, PS, plotcolor(rem(2*n+m,size(plotcolor,2))+1))
        ylabel('Amplitude')
        xlabel('Frequency(Hz)')
        title('Power Spectral of Dynamic FC between Seeds')
        legend('S1-S2','S1-S3','S1-S4','S2-S3','S2-S4','S3-S4')
        hold on
    end
end
saveas(figure(7),[figdir,'TC_dyn_FC_between_seeds.png']);
saveas(figure(8),[figdir,'PS_dyn_FC_between_seeds.png']);

%             % plot the time course of a state
%             figure(15+4*(i-1))
%             Time=ceil(W_width/2):step:(ceil(W_width/2)+step*(N_win-1));
%             plot(Time,indx)
%             ylim([0 8])
%             xlim([0 N_vol])
%             ylabel('States')
%             xlabel('Time')
%             title('Time Course of States')
%
%             % calculate percent occurence of each state
%             for j=1:N_stat
%                 N_occ(i,j)=sum(indx==j);
%                 perc_occ(i,j)=N_occ(j)/N_win;
%             end
%             save('/Users/zhenyang/Documents/microstate/Results/StatOcc.mat','N_occ','perc_occ');
%
%
%             % compute the transition matrix of Markov Chain
%             markovChain=indx;
%             Norder=1;
%             [ transitionMatrix columnStates ] = getTransitionMatrix(markovChain,Norder);
%             transitionMatrix=transitionMatrix'; %the row represent state at time t and each column the time t+1
%             figure(16+4*(i-1))
%             imagesc(transitionMatrix)
%             colorbar
%             ylabel('States at time t')
%             xlabel('States at time t+1')
%             title('Transition Matrix of Markov Chain')
%
%             [Eig_Vector,Eig_value] = eig(transitionMatrix);
%             stat_vector=Eig_Vector(:,1)';
%             figure(17+4*(i-1))
%             plot(stat_vector,'-or','Markersize', 10)
%             ylabel('Stationary Probability')
%             xlabel('States')
%             title('Steady-state Behavior')
%         end














