clear
clc
close all

% window width in TRs
W_width = 69;
% window sliding step in TRs
step=3;
% time
T=0:W_width-1;
% repetition time
TR=0.645;
% sampling rate
fs=1/TR;
% total time points
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
% extract the FC between each seed and ROIs
stat_FC=zeros(N_seed,N_ROI);
for i=1:N_seed
    stat_FC(i,:)=Corr(i,N_seed+1:N_ROI+N_seed);
end
% Fisher transformation
z=0.5*log((1+stat_FC)./(1-stat_FC));
% rearrange the regions by lobes
% z=z(i,[order of columns]);
% Plot the stationary FC matrix in Z
for i=1:N_seed
    figure(1)
    subplot(2,2,i)
    imagesc(z(i,:))
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

% Plot examplar dynamic FC matrix in r
for i=1:50:N_win
    figure(3)
    subplot(3,3,(i-1)/50+1)
    imagesc(Dyn_FC(:,:,i))
    colorbar('EastOutside')
    caxis([-1,1])
    set(gca,'YTickLabel','');
    ylabel('Correltion (r)')
    xlabel('ROIs')
    title('Dynamic FC')
end

% Plot examplar dynamic FC matrix in Z
for i=1:50:N_win
    figure(4)
    subplot(3,3,(i-1)/50+1)
    imagesc(z_Dyn_FC(:,:,i))
    colorbar('EastOutside')
    set(gca,'YTickLabel','');
    ylabel('Correltion (Z)')
    xlabel('ROIs')
    title('Dynamic FC')
end


% Examplar time course and power sepctral of dynamic FC

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
        title('Dynamic FC Time Course')
        hold on
        figure(i+8)
        t=-fs/2:fs/(length(PS)):fs/2-fs/(length(PS));
        plot(t, PS, plotcolor(rem(j,size(plotcolor,2))+1))
        ylabel('Amplitude')
        xlabel('Frequency(Hz)')
        title('Power Spectral of Dynamic FC Time Course')
        hold on
    end
end

% state identification

% exemplar selection
for k=1:N_win
    atemp=z_Dyn_FC(:,:,k)';
    btemp=var(atemp);
    for i=1:N_seed
        FC_var(k,i)=btemp(1,i);
    end
end
FC_var_Trans=FC_var';

% plot FC_var
figure(13)
for i=1:N_seed
    subplot(2,2,i)
    plot(FC_var_Trans(i,:));
    xlim([0 N_win]);
    ylabel('Variance')
    xlabel('Window Number')
    title('FC Variance')
end

for i=1:N_seed
    seed_FC_var=FC_var_Trans(i,:);
    [pks,locs] = findpeaks(seed_FC_var);
    N_exemplar=size(pks,2);
    for j=1:N_exemplar
        local_max(i,j)=pks(1,j);
        win_ID(i,j)=locs(1,j);
    end
    N_loc_max(i)=length(find(win_ID(i,:)>0));
end


% for each seed find and plot the state

for i=1:N_seed
    % select exemplar windows for clustering
    win_ID_final=win_ID(i,1:N_loc_max(i));
    for k=1:size(win_ID_final,2)
       loc_max_win(i,:,k)=z_Dyn_FC(1,:,win_ID_final(k));
    end
    
    % clear exemp_win indx ctrs all_win indx_final ctrs_final transitionMatrix Eig_Vector Eig_value 
    % k-mean clustering
    N_comp=7;
    exemp_win=reshape(loc_max_win(i,:,1:N_loc_max(i)),N_ROI,[])';
    [indx,ctrs]=kmeans(exemp_win,N_comp,'Distance','cityblock','Replicates',500);
    all_win=reshape(z_Dyn_FC(i,:,:),N_ROI,[])';
    [indx_final,ctrs_final]=kmeans(all_win,N_comp,'Distance','cityblock','emptyaction','singleton','start', ctrs);
    
    % calculate percent occurence of each state
    for j=1:N_comp
        N_occ(i,j)=sum(indx_final==j);
        perc_occ(i,j)=N_occ(j)/N_win;
    end
    save('/Users/zhenyang/Documents/microstate/Results/StatOcc.mat','N_occ','perc_occ');
    
    
    % plot the state
    figure(14+4*(i-1))
    imagesc(ctrs_final)
    colorbar('EastOutside')
    ylabel('State')
    xlabel('ROIs')
    title('Temporal Dynamics of FC')
    
    % plot the time course of a state
    figure(15+4*(i-1))
    Time=ceil(W_width/2):step:(ceil(W_width/2)+step*(N_win-1));
    plot(Time,indx_final)
    ylim([0 8])
    xlim([0 N_vol])
    ylabel('States')
    xlabel('Time')
    title('Time Course of States')
    
    % compute the transition matrix of Markov Chain
    markovChain=indx_final;
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

N_fig=29;
for i=1:N_fig
str1='fig';
filename=sprintf('%s%d',str1,i);
saveas(figure(i),['/Users/zhenyang/Documents/microstate/figures/69TR_step3/', filename,'.png']);
end



    %   % modularity analysis with Louvain algorithm
    %   for i=1:N_comp
    %   [Ci(i) Q(i)] = modularity_louvain_und_sign(ctrs1(i,:),'sta');
    %   end
    
    
    

  
  
   
    
    












