clear
clc

% N: window lenth in TRs
W_width = 22;
% T: time points
T = 0:W_width-1;
dr=60;
% sampling rate
TR=0.645
fs=1/TR;
N_vol=450;
N_win=N_vol-W_width+1;


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


% Z normalize the time series of each voxel
[AllVolume, VoxelSize, ImgFileList, Header1, nVolumn] =rest_to4d('/Users/zhenyang/Documents/microstate/data/0021018/swCovRegressed_4DVolume');
[nDim1 nDim2 nDim3 nDimTimePoints]=size(AllVolume);
BrainSize = [nDim1 nDim2 nDim3];

% Convert into 2D

MaskData=rest_loadmask(nDim1, nDim2, nDim3, '/Users/zhenyang/Documents/microstate/mask/mask_used/BrainMask_05_61x73x61.img');

MaskData =logical(MaskData);%Revise the mask to ensure that it contain only 0 and 1


AllVolume=reshape(AllVolume,[],nDimTimePoints)';

MaskDataOneDim=reshape(MaskData,1,[]);
MaskIndex = find(MaskDataOneDim);
 
AllVolume=AllVolume(:,MaskIndex);


% Z_norm the time series for each voxel
AllVolume = (AllVolume-repmat(mean(AllVolume),size(AllVolume,1),1))./repmat(std(AllVolume),size(AllVolume,1),1);   
AllVolume(isnan(AllVolume))=0;


% Convert into 4D
AllVolumeBrain = single(zeros(nDimTimePoints, nDim1*nDim2*nDim3));
AllVolumeBrain(:,MaskIndex) = AllVolume;
 
AllVolumeBrain=reshape(AllVolumeBrain',[nDim1, nDim2, nDim3, nDimTimePoints]);


% write 4D file as a nift file
rest_Write4DNIfTI(AllVolumeBrain,Header1,'/Users/zhenyang/Documents/microstate/data/0021018/norm_AllVolume.nii')

disp ('Time series of each voxel of Z score normalized.')

% Create seed and atlas masks

[Outdata1,VoxDim,Header2]=rest_readfile('/Users/zhenyang/Documents/microstate/mask/mask_used/mask_used_Resliced/Resliced_HarvardOxford-cort-maxprob-thr25-2mm_YCG.nii');
PreC_PCC=(Outdata1==300|Outdata1==301|Outdata1==310|Outdata1==311);
rest_WriteNiftiImage(PreC_PCC,Header2,'/Users/zhenyang/Documents/microstate/mask/mask_used/PreC_PCC.nii');

[Outdata2,VoxDim,Header3]=rest_readfile('/Users/zhenyang/Documents/microstate/mask/mask_used/mask_used_Resliced/Resliced_CC200ROI_tcorr05_2level_all.nii');

OutdataReduced=Outdata2;
OutdataReduced(find(PreC_PCC))=0;
rest_WriteNiftiImage(OutdataReduced,Header3,'/Users/zhenyang/Documents/microstate/mask/mask_used/reduced.nii');

for i=1:200
    if length(find(OutdataReduced==i))~=length(find(Outdata2==i))
       OutdataReduced(find(Outdata2==i))=0;
    else
       OutdataReduced=OutdataReduced;
    end
end
ROI_index=unique(OutdataReduced);
N_ROI=length(ROI_index);
rest_WriteNiftiImage(OutdataReduced,Header3,'/Users/zhenyang/Documents/microstate/mask/mask_used/final_reduced.nii');

disp ('Create a ROI mask on Crad-200 with Precuneus and PCC removed based on Havard-Oxford atlas.')

% extract Seed and ROI time series

def_seed={'ROI Center(mm)=(-2, -36, 35); Radius=3.00 mm.';'ROI Center(mm)=(-2, -47, 58); Radius=3.00 mm.';...
'ROI Center(mm)=(-2, -64, 45); Radius=3.00 mm.';'ROI Center(mm)=(-1, -78, 43); Radius=3.00 mm.'}

[seed_ROISignals] = y_ExtractROISignal('/Users/zhenyang/Documents/microstate/data/0021018/norm_AllVolume', def_seed,...
'/Users/zhenyang/Documents/microstate/data/time_series/seed_ROISignal',MaskData);

[atlas_ROISignals] = y_ExtractROISignal('/Users/zhenyang/Documents/microstate/data/0021018/norm_AllVolume', ...
{'/Users/zhenyang/Documents/microstate/mask/mask_used/final_reduced.nii'},...
'/Users/zhenyang/Documents/microstate/data/time_series/atlas_ROISignal',MaskData,1);

disp ('Extract the time series for the seeds and the selected ROIs.')


% For each seed, compute and plot the stationary FC matrix


seed_ROISignals = load('/Users/zhenyang/Documents/microstate/data/time_series/ROISignals_seed_ROISignal.mat');
TC1=seed_ROISignals.ROISignals;
N_seed=size(TC1,2);
%N_seed=1;
ROI_ROISignals=load('/Users/zhenyang/Documents/microstate/data/time_series/ROISignals_atlas_ROISignal.mat');
TC2=ROI_ROISignals.ROISignals;
N_ROI=size(TC2,2);
TC=[TC1,TC2];

% if calc_cor=1, generate full correlation matrix; if calc_cor=0,calc covariance and regularized precision matrix 
calc_cor = 1;

if (calc_cor == 1)
   Corr=corrcoef(TC);
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
else
   Cov=cov(TC);
   % comuput regularized covariance matrix
   % Fisher transformation
   z=0.5*log((1+stat_FC)./(1-stat_FC)); 
   % rearrange the regions by lobes
   %z=z(1,[order of columns]);
   % Plot the stationary FC matrix
   imagesc(z(i,:));
   ylabel('Correltion (Z)')
   xlabel('ROIs')
   title('Stationary FC')
end
   



% avrage across all subjects




% For each seed, compute and plot the dynamic FC matrix

% apply the sliding window to the time series
asize = size(TC);
N_win=asize(1)-W_width+1;
WC_final=WC(1,11:32)';

win=zeros(W_width,asize(2),N_win);
for k=1:N_win;
    for i=k:W_width+k-1
        for j=1:asize(2)
            win(i-k+1,j,k)=TC(i,j)*WC_final(i-k+1);
        end
    end
end
 
% if calc_cor=1, generate full correlation matrix; if calc_cor=0,calc covariance and regularized precision matrix 

if (calc_cor == 1)
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
   for i=1:50:401
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
   for i=1:50:401
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
       FC_TC_plot=[zeros(1,W_width/2),FC_TC];
       figure(i+4)
       plot(FC_TC_plot,plotcolor(rem(j,size(plotcolor,2))+1))
       xlim([0,N_vol])
       ylabel('Correltion (Z)')
       xlabel('Time(TR)')
       title('Dynamic FC Time Course')
       hold on
       figure(i+8)
       t=-fs/2:fs/(length(PS)):fs/2-fs/(length(PS))
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
  
  % find local maximum variance for each seed
  for i=1:N_seed
      seed_FC_var=FC_var_Trans(i,:);
      [pks,locs] = findpeaks(seed_FC_var);
      N_exemplar=size(pks,2);
      for j=1:N_exemplar
          local_max(i,j)=pks(1,j);
          win_ID(i,j)=locs(1,j);
      end
  end
  
   for i=1:N_seed
      N_loc_max(i)=length(find(win_ID(i,:)>0));
  end
  
  
  % seed 1
  % select exemplar windows for clustering
  win_ID_final=win_ID(1,1:N_loc_max(1));
  for k=1:size(win_ID_final,2)
      seed1_loc_max(1,:,k)=z_Dyn_FC(1,:,win_ID_final(k));
  end
  
  % k-mean clustering
  N_comp=7;
  seed1_exemplar=reshape(seed1_loc_max,N_ROI,[])';
  [indx,ctrs]=kmeans(seed1_exemplar,N_comp,'Distance','cityblock','Replicates',500);
  seed1_all_win=reshape(z_Dyn_FC(1,:,:),N_ROI,[])';
  [indx1,ctrs1]=kmeans(seed1_all_win,N_comp,'Distance','cityblock','start', ctrs);  

  % calculate percent occurence of each state
  for i=1:N_comp
  N_occ(i)=sum(indx1==i)
  perc_occ(i)=N_occ(i)/N_win
  end

  % plot the state
  figure(14)
  imagesc(ctrs1)
  colorbar('EastOutside')
  %set(gca,'YTickLabel','');
  ylabel('State')
  xlabel('ROIs')
  title('Temporal Dynamics of FC')
   
  % plot the time course of a state 
  figure(15)
  y=vertcat(zeros(W_width/2,1),indx1);
  plot(y)
  ylim([0 8])
  xlim([0 N_vol])
  
  [Eig_Vector,Eig_value] = eig(transitionMatrix);
  stat_vector=Eig_Vector(:,1)';
  figure(16)
  plot(stat_vector,'o')
  

%   % modularity analysis with Louvain algorithm
%   for i=1:N_comp
%   [Ci(i) Q(i)] = modularity_louvain_und_sign(ctrs1(i,:),'sta');
%   end

% compute the transition matrix of Markov Chain
markovChain=indx1;
Norder=1;
[ transitionMatrix columnStates ] = getTransitionMatrix(markovChain,Norder);
transitionMatrix=transitionMatrix'; %the row represent state at time t and each column the time t+1
figure(16)
imagesc(transitionMatrix)
colorbar
ylabel('States at time t')
xlabel('States at time t+1')
title('Transition Matrix of Markov Chain')


% seed 2
      win_ID_final=win_ID(2,1:N_loc_max(2));
      for k=1:size(win_ID_final,2)
          seed2_loc_max(1,:,k)=z_Dyn_FC(2,:,win_ID_final(k));
      end
      
      
   
  % seed 3
      win_ID_final=win_ID(3,1:N_loc_max(3));
      for k=1:size(win_ID_final,2)
          seed3_loc_max(1,:,k)=z_Dyn_FC(3,:,win_ID_final(k));
      end
  % seed 4
      win_ID_final=win_ID(4,1:N_loc_max(4));
      for k=1:size(win_ID_final,2)
          seed4_loc_max(1,:,k)=z_Dyn_FC(4,:,win_ID_final(k));
      end
      


else
   %generate the regularized covariance matrix
end

    
    
   
    
    
%     lamda
%    
%     GraphicalLassoPath(pop, lambdaList, approximate, verbose, penalDiag, tolThreshold, maxIter)
    











AM(abs(AM) < 1e-4) = 0;
AM2 = sign(AM);
drawGraph(abs(AM2),'labels',wordlist);


ts1=W1;


lambda=5; % arbitrary choice of regularisation!
oc=cov(ts1); % raw covariance
ic=-L1precisionBCD(oc/mean(diag(oc)),lambda/1000); % get regularised negative inverse covariance
r=(ic ./ repmat(sqrt(abs(diag(ic))),1,Nnodes)) ./ repmat(sqrt(abs(diag(ic)))',Nnodes,1); % use diagonal to get normalised coefficients
r=r+eye(Nnodes); % remove diagonal 






% applying the sliding window to the time course


% % plot the Gaussian window
% figure
% area(k,w2,'FaceColor', [0 .4 .6])
% xlim([0 N-1])
% set(gca,'XTick', [0 : 1/8 : 1]*(N-1))
% % set(gca,'XTickLabel','0| | | | | | | |N-1')
% grid on
% ylabel('amplitude')
% xlabel('Time')
% title('Window function (Gauss, \sigma = )')

% Data1(Data2==99)=0;
% TC_seed = load ('data/Seed_time_series/FunRaw_ROISignals.mat')
%concatenate files
%ROISignals = [seed_ROISignals(:,i),atlas_ROISignals];

% rearrange the regions by lobes
   %z=z([order of rows],[order of columns]);


% % evaluate the freq response of WR
% H1 = abs(fft(RW));
% H1 = fftshift(H1);
% H1 = H1/max(H1);
% H1 = 20*log10(H1);
% H1 = max(-dr,H1);
% h1=-fs/2:fs/(length(H1)):fs/2-fs/(length(H1))
% figure(1); subplot(1,2,2); plot(h1,H1)
% xlim([-0.25,0.25])
% 
% 
% % evaluate the freq response of WG
% H2 = abs(fft(GW));
% H2 = fftshift(H2);
% H2 = H2/max(H2);
% H2 = 20*log10(H2);
% H2 = max(-dr,H2);
% h2=-fs/2:fs/(length(H2)):fs/2-fs/(length(H2));
% figure(2); subplot(1,2,2);plot(h2,H2)
% xlim([-0.25,0.25])
% 
% % evaluate the freq response of WC
% H3= abs(fftshift(fft(WC)));
% H3= H3/max(H3);
% H3= 20*log10(H3);
% H3= max(-dr,H3);
% h3=-fs/2:fs/(length(H3)):fs/2-fs/(length(H3));
% figure(3); subplot(1,2,2);plot(h3,H3)
% xlim([-0.25,0.25])
% 
% % creat the hanning window and hamming window and evaluate their frequency
% % responses
% HanW = hann(W_width+1)
% figure(4); subplot(1,2,1); plot(HanW)
% xlim([0 length(HanW)-1])
% 
% H4= abs(fftshift(fft(HanW)));
% H4= H4/max(H4);
% H4= 20*log10(H4);
% H4= max(-dr,H4);
% h4=-fs/2:fs/(length(H4)):fs/2-fs/(length(H4));
% figure(4); subplot(1,2,2);plot(h4,H4)
% xlim([-0.25,0.25])
% 
% HammingW = hamming(W_width+1)
% figure(5); subplot(1,2,1); plot(HammingW)
% xlim([0 length(HammingW)-1])
% 
% H5= abs(fftshift(fft(HammingW)));
% H5= H5/max(H5);
% H5= 20*log10(H5);
% H5= max(-dr,H5);
% h5=-fs/2:fs/(length(H5)):fs/2-fs/(length(H5));
% figure(5); subplot(1,2,2);plot(h5,H5)
% xlim([-0.25,0.25])
% 
% 
% 
% plot(y, 'color', [0.4 0.5 0.6]) RGB value between 0~1
% 
% 
