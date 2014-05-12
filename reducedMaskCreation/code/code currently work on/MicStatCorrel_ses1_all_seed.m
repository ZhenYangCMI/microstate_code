clear
clc
close all


subList={'0021002', '0021006', '0021018', '0021024', '1427581', '1793622',...
    '1961098', '2475376', '2799329', '2842950', '3201815', '3313349'...
    '3315657', '3795193', '3808535', '3893245', '4176156', '4288245',...
    '7055197', '8574662', '8735778', '9630905'};
%SesList={'session1','session2'};
%TRList={'645','2500'};

% SubList={'0021002', '0021006'};
sesList={'session2'};
TRList={'645'};

numSub=size(subList,2);
numSes=size(sesList,2);
numTR=size(TRList,2);
numSeed=4;
numROI=156;


ROIIndx=load('/home/data/Projects/microstate/analysis_matlab/doc/Crad179_ROI_index.txt');
% rearrange the regions by lobes, after rearrange:lobe1=F (75ROIs, col 1:75); lobe2=T (31ROIs, col 76:106); lobe3=p(20ROIs,col 107:126); lobe4=o(21ROIs,col 127:147);
% lobe5=subcortical(9ROIs, col 148:156)
orderByLobe=xlsread('/home/data/Projects/microstate/analysis_matlab/doc/Crad_179_lob_hier_label_final.xls','Crad179_ROI_index','a2:c157');
orderByLobe=orderByLobe(:,3)';
analyDir=['/home/data/Projects/microstate/analysis_matlab/'];
figDir=['/home/data/Projects/microstate/analysis_matlab/fig/'];
maskDir=['/home/data/Projects/microstate/analysis_matlab/mask/'];

for i=1:numTR
    TR=TRList{i};
    TRDir=['/home/data/Projects/microstate/analysis_matlab/data/',TR];
    if ~exist(['/home/data/Projects/microstate/analysis_matlab/results/',TR], 'dir')
        mkdir('/home/data/Projects/microstate/analysis_matlab/results', TR)
    end
    resultTRDir=['/home/data/Projects/microstate/analysis_matlab/results/',TR];
    for j=1:numSes
        ses=sesList{j};
        sesDir=[TRDir,'/',ses];
        if ~exist([resultTRDir,'/',ses], 'dir')
            mkdir(resultTRDir, ses)
        end
        resultSesDir=[resultTRDir,'/',ses];
        
        if(~exist([resultSesDir,'/zFCMean.mat'],'file')
            for k=1:numSub
                subDir=[sesDir,'/', char(subList{k})];
                disp (['Working on sub ', char(subList{k}),' ......'])
                
                seedROISignals = load([subDir,'/ROISignals_seed_ROISignal.mat']);
                TC1=seedROISignals.ROISignals;
                numSeed=size(TC1,2);
                
                ROIROISignals=load([subDir,'/ROISignals_atlas_ROISignal.mat']);
                TC2=ROIROISignals.ROISignals;
                numROI=size(TC2,2);
                
                % concatenate the time series of seeds and ROIs
                TC=[TC1,TC2];
                
                % Stationary FC estimation
                
                % Compute the full correlation matrix
                stateCor=zeros(size(TC,2),size(TC,2),N_sub);
                stateCor(:,:,k)=corrcoef(TC);
                
                % Extract the stationary FC between all pairs of seeds
                for n=1:numSeed
                    for m=1:numSeed
                        FCBtwSeed(n,m,k)=stateCor(n,m,k);
                        zFCBetwSeed(n,m,k)=0.5*log((1+FCBtwSeed(n,m,k))./(1-FCBtwSeed(n,m,k)));
                    end
                end
                
                % Extract the stationary FC between seeds and ROIs
                FC=zeros(numSeed,numROI,numSub);
                for n=1:numSeed
                    for m=numSeed+1:numROI+numSeed
                        FC(n,m-numSeed,k)=stateCor(n,m,k);
                        zFC(n,m-numSeed,k)=0.5*log((1+FC(n,m-numSeed,k))./(1-FC(n,m-numSeed,k)));
                    end
                end
                
            end
            
            % average across all subjects
            temp1=reshape(zFCBetwSeed,[],k)';
            temp1Mean=mean(temp1);
            zFCBtwSeedMean=reshape(temp1Mean,numSeed,numSeed);
            save([ses_result_dir,'/zFCBtwSeedMean.mat'],'zFCBtwSeedMean');
            
            temp2=reshape(zFC,[],k)';
            temp2Mean=mean(temp2);
            zFCMean=reshape(temp2Mean,numSeed,numROI);
            save([ses_result_dir,'/zFCMean.mat'],'zFCMean');
            
            % plot the stat_FC between seeds
            figure(1)
            imagesc(zFCBtwSeedMean)
            colorbar('EastOutside')
            set(gca,'XTick',1:4,'YTick',1:4);
            ylabel('Seeds')
            xlabel('Seeds')
            %xlim([-1 1]);
            title('Stat FC between Seeds (Z)')
            saveas(figure(1),[figdir, 'stat_FC_between_seeds.png']);
            
            % rearrange stat_FC ROIs by lobe
            
            zFCMeanByLobe=zFCMean(:,[col_order]);
            
            % Plot the stat_FC(Z) between seeds and ROIs
            figure(2)
            imagesc(zFCMeanByLobe)
            colorbar('EastOutside')
            caxis([-0.5,0.5])
            set(gca,'YTickLabel','');
            ylabel('Seeds')
            xlabel('ROIs')
            title('Stat FC between Seeds and ROIs (Z)')
            saveas(figure(2),[figdir, 'stat_FC_seeds_and_ROIs.png']);
            
            disp('Stat_FC estimation done!')
        end
        
        
        
        
        % Dynamic FC estimation
        
        for k=1:numSub
            subDir=[sesdir,'/', char(SubList{k})];
            disp (['Working on sub ', char(SubList{k}),' ......'])
            
            seedROISignals = load([subDir,'/ROISignals_seed_ROISignal.mat']);
            TC1=seedROISignals.ROISignals;
            numSeed=size(TC1,2);
            
            ROIROISignals=load([subDir,'/ROISignals_atlas_ROISignal.mat']);
            TC2=ROIROISignals.ROISignals;
            numROI=size(TC2,2);
            
            % concatenate the time series of seeds and ROIs
            TC=[TC1,TC2];
            
            
            % apply the sliding window to the time series
            asize = size(TC);
            [finalWin]=createWin(winWidth);
            finalWin=finalWin';
            % win(W_width,asize(2),N_win_tot);
            for q=1+numWin*(k-1):numWin*k;
                for n=((q-1)-(k-1)*numWin)*step+1:winWidth+((q-1)-(k-1)*numWin)*step
                    for m=1:asize(2)
                        win(n-((q-1)-(k-1)*numWin)*step,m,q)=TC(n,m)*finalWin((n-((q-1)-(k-1)*numWin)*step),1);
                    end
                end
            end
        end
        save([resultSesDir,'/win_all_sub.mat'],'win');
        
        disp ('Window applying done!')
        
        % generate the full correlation matrix for each win and extrac the dyn FC
        
        % Dyn_FC(N_seed,N_ROI,numWin_tot);
        % Dyn_FC_seed(N_seed,N_seed,numWin_tot);
        numWinAllSub=numWin*numSub;
        for q=1:numWinAllSub
            corrWin(:,:,q)=corrcoef(win(:,:,q));
            for n=1:numSeed
                for m=1:numSeed
                    dynFCBtwSeed(n,m,q)=corrWin(n,m,q);
                    zDynFCBtwSeed(n,m,q)=0.5*log((1+dynFCBtwSeed(n,m,q))./(1-dynFCBtwSeed(n,m,q)));
                end
                for m=numSeed+1:numROI+numSeed
                    dynFC(n,m-numSeed,q)=dynCorr(n,m,q);
                    zDynFC(n,m-numSeed,q)=0.5*log((1+dynFC(n,m-numSeed,q))./(1-dynFC(n,m-numSeed,q)));
                end
            end
        end
        save([ses_result_dir,'/Dyn_FC.mat'],'Dyn_FC');
        save([ses_result_dir,'/z_Dyn_FC.mat'],'z_Dyn_FC');
        disp ('Full correlation of each window computed!')
    end
end

loadfile=load([resultSesDir,'/z_Dyn_FC.mat'])
zDynFC=loadfile.z_Dyn_FC;


N_stat=7;


win_all_sub1=reshape(zDynFC(1,:,:),numROI,[])';
win_all_sub2=reshape(zDynFC(2,:,:),numROI,[])';
win_all_sub3=reshape(zDynFC(3,:,:),numROI,[])';
win_all_sub4=reshape(zDynFC(4,:,:),numROI,[])';
win_all_sub_seed=vertcat(win_all_sub1, win_all_sub2, win_all_sub3, win_all_sub4);
disp ('Win of all seeds concatenated.')
for N_clust=2:20
    clear indx_stat ctrs_stat sumd D
    opts=statset('MaxIter',500);
    [indx_stat,ctrs_stat,sumd,D]=kmeans(win_all_sub_seed,N_stat,'Distance','correlation','emptyaction','singleton','Replicates',200,'options',opts);
    % compute the cluster validity index
    CVI(N_clust)=sum(sumd)./(sum(sum(D))-sum(sumd));
    %         index_stat(:,n)=indx_stat;
    %         final_ctrs(:,:,n)=ctrs_stat;
    %         final_ctrs_transp(:,:,n)=ctrs_stat';
    disp ('k-means clustering with',N_clust,' cluster done.')
end
%         final_ctrs_cat=reshape(final_ctrs_transp, N_ROI, []);
final_ctrs=ctrs_stat;
final_ctrs_transp=ctrs_stat';
% compare the similarity between states
[stat_sim,p_val]=corrcoef(final_ctrs_transp);
save([analyDir,'results/states_all_sub_all_seeds.txt'],'ctrs_stat','-ascii', '-double', '-tabs')
save([analyDir,'results/state_sim.txt'],'stat_sim','-ascii', '-double', '-tabs')
save([analyDir,'results/p_val.txt'],'p_val','-ascii', '-double', '-tabs')
save([analyDir,'results/indx_all_sub_all_seeds.txt'],'indx_stat','-ascii', '-double', '-tabs')
save([analyDir,'results/sumd_all_sub_all_seeds.txt'],'sumd','-ascii', '-double', '-tabs')
save([analyDir,'results/btD_all_sub_all_seeds.txt'],'D','-ascii', '-double', '-tabs')

% plot the correlation of states between seeds
figure(9)
imagesc(stat_sim)
colorbar('EastOutside')
caxis([-1,1])
set(gca,'XTick',1:3:numSeed*N_stat,'YTick', 1:3:numSeed*N_stat);
ylabel('States of all seeds')
xlabel('States of all seeds')
title('Correlations of states of all seeds')
saveas(figure(9),[figDir,'corr_stat_all_seeds.png']);

%     figure(10)
%     subplot(3,2,1)
%     imagesc(stat_sim(N_stat+1:N_stat*2,1:N_stat))
%     colorbar('EastOutside')
%     caxis([-1,1])
%     set(gca,'XTick',1:N_stat,'YTick', 1:N_stat);
%     ylabel('Seed 2 States')
%     xlabel('Seed 1 States')
%     title('State Similarity')
%
%     subplot(3,2,2)
%     imagesc(stat_sim(N_stat*2+1:N_stat*3,1:N_stat))
%     colorbar('EastOutside')
%     caxis([-1,1])
%     set(gca,'XTick',1:N_stat,'YTick', 1:N_stat);
%     ylabel('Seed 3 States')
%     xlabel('Seed 1 States')
%     title('State Similarity')
%
%     subplot(3,2,3)
%     imagesc(stat_sim(N_stat*3+1:N_stat*4,1:N_stat))
%     colorbar('EastOutside')
%     caxis([-1,1])
%     set(gca,'XTick',1:N_stat,'YTick', 1:N_stat);
%     ylabel('Seed 4 States')
%     xlabel('Seed 1 States')
%     title('State Similarity')
%
%     subplot(3,2,4)
%     imagesc(stat_sim(N_stat*2+1:N_stat*3,N_stat+1:N_stat*2))
%     colorbar('EastOutside')
%     caxis([-1,1])
%     set(gca,'XTick',1:N_stat,'YTick', 1:N_stat);
%     ylabel('Seed 3 States')
%     xlabel('Seed 2 States')
%     title('State Similarity')
%
%     subplot(3,2,5)
%     imagesc(stat_sim(N_stat*3+1:N_stat*4,N_stat+1:N_stat*2))
%     colorbar('EastOutside')
%     caxis([-1,1])
%     set(gca,'XTick',1:N_stat,'YTick', 1:N_stat);
%     ylabel('Seed 4 States')
%     xlabel('Seed 2 States')
%     title('State Similarity')
%
%     subplot(3,2,6)
%     imagesc(stat_sim(N_stat*3+1:N_stat*4,N_stat*2+1:N_stat*3))
%     colorbar('EastOutside')
%     caxis([-1,1])
%     set(gca,'XTick',1:N_stat,'YTick', 1:N_stat);
%     ylabel('Seed 4 States')
%     xlabel('Seed 3 States')
%     title('State Similarity')
%     saveas(figure(10),[figdir,'corr_stat_seed_pair.png']);

% plot the state
figure (11)
for n=1:numSeed
    subplot(2,2,n)
    final_ctrs_bylobe=final_ctrs(:,[orderByLobe]);
    imagesc(final_ctrs_bylobe(:,:))
    colorbar('EastOutside')
    ylabel('States')
    xlabel('ROIs')
    title(['Seed',num2str(n),'States'])
end
saveas(figure(11),[figDir,'states.png']);

[Outdata,VoxDim,Header]=rest_readfile([maskDir,'final_reduced.nii']);
[nDim1 nDim2 nDim3]=size(Outdata);
temp=unique(Outdata);
ROIIndx=temp(find(temp~=0));
% for n=1:N_seed
for t=1:N_stat
    Statemap=Outdata;
    state_ctrs=final_ctrs(t,:);
    for m=1:numROI
        Statemap(find(Outdata==ROIIndx(m)))=state_ctrs(m);
    end
    
    Header.pinfo = [1;0;0];
    Header.dt    =[16,0];
    rest_WriteNiftiImage(Statemap,Header,[figDir,'state', num2str(t),'.nii']);
end
%end
end
disp([ses,' done!'])
end

