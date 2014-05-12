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
% Gaussian Standard Deviation
Gs=3;

GW=exp(-(T-Gm).^2/(2*Gs.^2));

% figure(2); subplot(1,2,1); plot(GW)
% xlim([0 length(GW)-1])

% convolve rectangular and gaussian
WC=conv(RW,GW);
% figure(3); subplot(1,2,1); plot(WC)
% xlim([0 length(WC)-1])


% For each seed, compute and plot the stationary FC matrix

seed_ROISignals = load('/Users/zhenyang/Documents/microstate/data/time_series/ROISignals_seed_ROISignal.mat');
TC1=seed_ROISignals.ROISignals;
N_seed=size(TC1,2);
%N_seed=1;
ROI_ROISignals=load('/Users/zhenyang/Documents/microstate/data/time_series/ROISignals_atlas_ROISignal.mat');
TC2=ROI_ROISignals.ROISignals;
N_ROI=size(TC2,2);
TC=[TC1,TC2];


Corr=corrcoef(TC);
stat_FC=zeros(N_seed,N_ROI);
for i=1:N_seed
    stat_FC(i,:)=Corr(i,N_seed+1:N_ROI+N_seed);
end
% Fisher transformation
z_stat_FC=0.5*log((1+stat_FC)./(1-stat_FC));
% rearrange the regions by lobes
% z_stat_FC_bylobe=z_stat_FC(i,[order of columns]);
% Plot the stationary FC matrix in Z
for i=1:N_seed
    figure(1)
    subplot(2,2,i)
    imagesc(z_stat_FC(i,:))
    colorbar('EastOutside')
    caxis([-1,1])
    set(gca,'YTickLabel','');
    ylabel('Correltion (Z)')
    xlabel('ROIs')
    title('Stationary FC')
    % Plot the stationary FC matrix in r
    figure(2)
    subplot(2,2,i)
    imagesc(stat_FC(i,:))
    colorbar('EastOutside')
    caxis([-1,1])
    set(gca,'YTickLabel','');
    ylabel('Correltion (r)')
    xlabel('ROIs')
    title('Stationary FC')
end
 
% extract the FC between all pairs of seeds
for i=1:N_seed
    stat_seed_FC(i,:)=Corr(i,1:N_seed);
end
% Fisher transformation
z_seed=0.5*log((1+stat_seed_FC)./(1-stat_seed_FC));
for i=1:N_seed
    figure(30)
    subplot(2,2,i)
    imagesc(z_seed(i,:))
    colorbar('EastOutside')
    caxis([-1,1])
    set(gca,'YTickLabel','');
    ylabel('Correltion (Z)')
    xlabel('Seeds')
    title('Stationary Seed FC')
end
    
% avrage across all subjects




% For each seed, compute and plot the dynamic FC matrix

% apply the sliding window to the time series
asize = size(TC);
WC_final=WC(1,ceil(W_width/2):(W_width+ceil(W_width/2)-1))';
win=zeros(W_width,asize(2),N_win);
for k=1:N_win;
    for i=(k-1)*step+1:W_width+(k-1)*step 
        for j=1:asize(2)
            win(i-(k-1)*step,j,k)=TC(i,j)*WC_final(i-(k-1)*step);
        end
    end
end



% generate the full correlation matrix

Dyn_FC=zeros(N_seed,N_ROI,N_win);
for k=1:N_win
    Corr=corrcoef(win(:,:,k));
    for i=1:N_seed
        for j=N_seed+1:N_ROI+N_seed
            Dyn_FC(i,j-N_seed,k)=Corr(i,j);
            z_Dyn_FC(i,j-N_seed,k)=0.5*log((1+Dyn_FC(i,j-N_seed,k))./(1-Dyn_FC(i,j-N_seed,k)));
        end
    end
end

for k=1:N_win
    Corr=corrcoef(win(:,:,k));
    for i=1:N_seed
        for j=1:N_seed
            Dyn_FC_seed(i,j,k)=Corr(i,j);
            z_Dyn_FC_seed(i,j,k)=0.5*log((1+Dyn_FC_seed(i,j,k))./(1-Dyn_FC_seed(i,j,k)));
        end
    end
end


% Plot examplar dynamic FC between seeds and ROIs in r
for i=1:50:N_win
    figure(3)
    subplot(3,3,(i-1)/50+1)
    imagesc(Dyn_FC(:,:,i))
    colorbar('EastOutside')
    caxis([-1,1])
    set(gca,'YTickLabel','');
    ylabel('Correltion (r)')
    xlabel('ROIs')
    title('Dynamic FC between Seeds and ROIs')
end

% Plot examplar dynamic FC between seeds and ROIs in Z
for i=1:50:N_win
    figure(4)
    subplot(3,3,(i-1)/50+1)
    imagesc(z_Dyn_FC(:,:,i))
    colorbar('EastOutside')
    set(gca,'YTickLabel','');
    ylabel('Correltion (Z)')
    xlabel('ROIs')
    title('Dynamic FC between Seeds and ROIs')
end


% Examplar time course and power sepctral of dynamic FC between seeds and ROIs

plotcolor=['r' 'g' 'b' 'c' 'm' 'y' 'k'];

for i=1:N_seed
    for j=1:50:N_ROI;
        FC_TC=reshape(z_Dyn_FC(i,j,:),1,[]);
        PS= abs(fftshift(fft(FC_TC)));
        figure(i+4)
        Time=ceil(W_width/2):step:(ceil(W_width/2)+step*(N_win-1));
        plot(Time,FC_TC,plotcolor(rem(j,size(plotcolor,2))+1))
        xlim([0,N_vol])
        ylabel('Correltion (Z)')
        xlabel('Time(TR)')
        title('Time Course of Dynamic FC between Seeds and ROIs')
        hold on
        figure(i+8)
        t=-fs/2:fs/(length(PS)):fs/2-fs/(length(PS));
        plot(t, PS, plotcolor(rem(j,size(plotcolor,2))+1))
        ylabel('Amplitude')
        xlabel('Frequency(Hz)')
        title('Power Spectral of Dynamic FC between Seeds and ROIs')
        hold on
    end
end
%
% Plot examplar dynamic FC among seeds in r
for i=1:50:N_win
    figure(13)
    subplot(3,3,(i-1)/50+1)
    imagesc(Dyn_FC_seed(:,:,i))
    colorbar('EastOutside')
    caxis([-1,1])
    set(gca,'YTickLabel','');
    ylabel('Correltion (r)')
    xlabel('Seeds')
    title('Dynamic FC among Seeds')
end

% Plot examplar dynamic FC matrix among seeds in Z
for i=1:50:N_win
    figure(14)
    subplot(3,3,(i-1)/50+1)
    imagesc(z_Dyn_FC_seed(:,:,i))
    colorbar('EastOutside')
    set(gca,'YTickLabel','');
    ylabel('Correltion (Z)')
    xlabel('Seeds')
    title('Dynamic FC among Seeds')
end


% Examplar time course and power sepctral of dynamic FC among seeds

plotcolor=['r' 'g' 'b' 'c' 'm' 'y' 'k'];

for i=1:N_seed
    for j=1:N_seed;
        FC_seed_TC=reshape(z_Dyn_FC_seed(i,j,:),1,[]);
        PS= abs(fftshift(fft(FC_seed_TC)));
        figure(i+14)
        Time=ceil(W_width/2):step:(ceil(W_width/2)+step*(N_win-1));
        plot(Time,FC_seed_TC,plotcolor(rem(j,size(plotcolor,2))+1))
        xlim([0,N_vol])
        ylabel('Correltion (Z)')
        xlabel('Time(TR)')
        title('Time Course of Dynamic FC among Seeds')
        hold on
        figure(i+18)
        t=-fs/2:fs/(length(PS)):fs/2-fs/(length(PS));
        plot(t, PS, plotcolor(rem(j,size(plotcolor,2))+1))
        ylabel('Amplitude')
        xlabel('Frequency(Hz)')
        title('Power Spectral of Dynamic FC among Seeds')
        hold on
    end
end

% for each seed find and plot the state

% delete the state map generated before running the following loop
% delete 
for i=1:N_seed
        
    % k-mean clustering
    % for N_comp=2:10
    N_comp=7;
    all_win=reshape(z_Dyn_FC(i,:,:),N_ROI,[])';
    [indx,ctrs,sumd,D]=kmeans(all_win,N_comp,'Distance','cityblock','Replicates',500);
   
       
    for j=1:N_comp
        for k=1:N_ROI
            final_ctrs(j,k,i)=ctrs(j,k);
        end
    end
    
    final_ctrs_cat=[final_ctrs(:,:,1)',final_ctrs(:,:,2)', final_ctrs(:,:,3)',final_ctrs(:,:,4)'];
    [state_sim,p]=corrcoef(final_ctrs_cat); % compare the similarity among states
    
    % compute the cluster validity index
    CVI=sum(sumd)./(sum(sum(D))-sum(sumd));
    %plot(CVI)
    
    % calculate percent occurence of each state
    for j=1:N_comp
        N_occ(i,j)=sum(indx==j);
        perc_occ(i,j)=N_occ(j)/N_win;
    end
    save('/Users/zhenyang/Documents/microstate/Results/StatOcc.mat','N_occ','perc_occ');
    
    N_state=7;
    
    % plot the state
    figure(14+4*(i-1))
    imagesc(ctrs)
    colorbar('EastOutside')
    ylabel('State')
    xlabel('ROIs')
    title('Temporal Dynamics of FC')
    
    [Outdata,VoxDim,Header]=rest_readfile('/Users/zhenyang/Documents/microstate/mask/final_ROIs_Crad179/final_reduced.nii');
    [nDim1 nDim2 nDim3]=size(Outdata);
    temp=unique(Outdata);
    ROI_index=temp(find(temp~=0));
    for j=1:N_state
        Statemap=Outdata;
        state_ctrs=ctrs(j,:);
        for k=1:N_ROI
            Statemap(find(Outdata==ROI_index(k)))=state_ctrs(k);
        end

        Header.pinfo = [1;0;0];
        Header.dt    =[16,0];
        rest_WriteNiftiImage(Statemap,Header,['/Users/zhenyang/Documents/microstate/figures/state_maps/state_', num2str(j),]);
    end
    
    % plot the time course of a state
    figure(15+4*(i-1))
    Time=ceil(W_width/2):step:(ceil(W_width/2)+step*(N_win-1));
    plot(Time,indx)
    ylim([0 8])
    xlim([0 N_vol])
    ylabel('States')
    xlabel('Time')
    title('Time Course of States')
    
    % compute the transition matrix of Markov Chain
    markovChain=indx;
    Norder=1;
    [ transitionMatrix columnStates ] = getTransitionMatrix(markovChain,Norder);
    transitionMatrix=transitionMatrix'; %the row represent state at time t and each column the time t+1
    figure(16+4*(i-1))
    imagesc(transitionMatrix)
    colorbar
    ylabel('States at time t')
    xlabel('States at time t+1')
    title('Transition Matrix of Markov Chain')
    
    [Eig_Vector,Eig_value] = eig(transitionMatrix);
    stat_vector=Eig_Vector(:,1)';
    figure(17+4*(i-1))
    plot(stat_vector,'-or','Markersize', 10)
    ylabel('Stationary Probability')
    xlabel('States')
    title('Steady-state Behavior')
end





% N_fig=29;
% for i=1:N_fig
% str1='fig';
% str2='/Users/zhenyang/Documents/microstate/figures/';
% str3='TR_step';
% str4='/';
% filepath=sprintf('%s%d%s%d%s',str2,W_width,str3,step,str4)
% filename=sprintf('%s%d',str1,i);
% saveas(figure(i),[filepath, filename,'.png']);
% end



    %   % modularity analysis with Louvain algorithm
    %   for i=1:N_comp
    %   [Ci(i) Q(i)] = modularity_louvain_und_sign(ctrs1(i,:),'sta');
    %   end
    

  
  
   
    
    












