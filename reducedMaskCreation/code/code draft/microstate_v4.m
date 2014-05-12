clear
clc

% N: window lenth in TRs
W_width = 22;
% T: time points
T = 0:W_width-1;
dr=60;
% sampling rate
fs=0.5;
N_vol=450;
N_win=N_vol-W_width+1;


% creat the rectangular window
RW = ones(1,W_width);
figure(1); subplot(1,2,1); plot(T,RW)
xlim([0 length(RW)-1])

% Calculate gaussian 
% Gaussian Mean
Gm=(W_width-1)/2;
% Gaussian Standard Deviation
Gs=3;

GW=exp(-(T-Gm).^2/(2*Gs.^2));

figure(2); subplot(1,2,1); plot(GW)
xlim([0 length(GW)-1])

% convolve rectangular and gaussian
WC=conv(RW,GW);
figure(3); subplot(1,2,1); plot(WC)
xlim([0 length(WC)-1])

% evaluate the freq response of WR
H1 = abs(fft(RW));
H1 = fftshift(H1);
H1 = H1/max(H1);
H1 = 20*log10(H1);
H1 = max(-dr,H1);
h1=-fs/2:fs/(length(H1)):fs/2-fs/(length(H1))
figure(1); subplot(1,2,2); plot(h1,H1)
xlim([-0.25,0.25])


% evaluate the freq response of WG
H2 = abs(fft(GW));
H2 = fftshift(H2);
H2 = H2/max(H2);
H2 = 20*log10(H2);
H2 = max(-dr,H2);
h2=-fs/2:fs/(length(H2)):fs/2-fs/(length(H2));
figure(2); subplot(1,2,2);plot(h2,H2)
xlim([-0.25,0.25])

% evaluate the freq response of WC
H3= abs(fftshift(fft(WC)));
H3= H3/max(H3);
H3= 20*log10(H3);
H3= max(-dr,H3);
h3=-fs/2:fs/(length(H3)):fs/2-fs/(length(H3));
figure(3); subplot(1,2,2);plot(h3,H3)
xlim([-0.25,0.25])

% creat the hanning window and hamming window and evaluate their frequency
% responses
HanW = hann(W_width+1)
figure(4); subplot(1,2,1); plot(HanW)
xlim([0 length(HanW)-1])

H4= abs(fftshift(fft(HanW)));
H4= H4/max(H4);
H4= 20*log10(H4);
H4= max(-dr,H4);
h4=-fs/2:fs/(length(H4)):fs/2-fs/(length(H4));
figure(4); subplot(1,2,2);plot(h4,H4)
xlim([-0.25,0.25])

HammingW = hamming(W_width+1)
figure(5); subplot(1,2,1); plot(HammingW)
xlim([0 length(HammingW)-1])

H5= abs(fftshift(fft(HammingW)));
H5= H5/max(H5);
H5= 20*log10(H5);
H5= max(-dr,H5);
h5=-fs/2:fs/(length(H5)):fs/2-fs/(length(H5));
figure(5); subplot(1,2,2);plot(h5,H5)
xlim([-0.25,0.25])

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
ROI_TC=TC2;
N_ROI=size(ROI_TC,2);

% if calc_cor=1, generate full correlation matrix; if calc_cor=0,calc covariance and regularized precision matrix 
calc_cor = 1;

for i=1:N_seed
   if (calc_cor == 1)
   seed_TC=TC1(:,i);
   for j=1:N_ROI
    temp=corrcoef(seed_TC,ROI_TC(:,j));
    stat_FC(i,j)=temp(1,2);
   end
   % Fisher transformation
   z=0.5*log((1+stat_FC)./(1-stat_FC)); 
   % rearrange the regions by lobes
   % z=z(i,[order of columns]);
   % Plot the stationary FC matrix in Z
   figure(1)
   subplot(2,2,i)
   imagesc(z(i,:))
   colorbar('EastOutside')
   set(gca,'YTickLabel','');
   ylabel('Correltion (Z)')
   xlabel('ROIs')
   title('Stationary FC')
   % Plot the stationary FC matrix in r
   figure(2)
   subplot(2,2,i)
   imagesc(stat_FC(i,:))
   colorbar('EastOutside')
   set(gca,'YTickLabel','');
   ylabel('Correltion (r)')
   xlabel('ROIs')
   title('Stationary FC')
   else
   TC=[seed_TC,ROI_TC];
   stat_FC=cov(TC);
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
end    




% avrage across all subjects




% For each seed, compute and plot the dynamic FC matrix



for i=1:N_seed
TC=[TC1,TC2]
TC2
    
% segment time course into windows
asize = size(TC1);
bsize = size(TC2);
WC=WC(1,11:32)';
WC1=repmat(WC,W_width,asize(2));
WC2=repmat(WC,W_width,bsize(2));
N_win=asize(1)-W_width+1;
for i = 1:N_win
    % apply the sliding window to the time series
    seed_win=zeros(W_width,
    seed_win=WC1.*TC1(i:W_width-1+i,:);
    ROI_win_1=WC2.*TC2(i:W_width-1+i,:)
    
    
    % creat the covraince matrix
    Cov=cov(TC_win);
    s=['cov_' int2str(i) '=' 'Cov'];
    eval(s)    
    str1='cov_';
    atemp=eval(sprintf('%s%d',str1,i));
    
    % create regualrized invese matrix
    
%     lamda
%    
%     GraphicalLassoPath(pop, lambdaList, approximate, verbose, penalDiag, tolThreshold, maxIter)
    
    % Fisher transform
    z=0.5*log((1+atemp)./(1-atemp));
    s=['z' int2str(i) '=' 'z'];
    eval(s)    
    str1='z';
    btemp=eval(sprintf('%s%d',str1,i))
    btemp(:,:,i)= btemp;
           
end


% Examplar time course of FC dynamics
ROI1=1
ROI2=2
ROI3=3
ROI4=4
Time=W_width/2:2:450
pair1=btemp(ROI1,ROI2,:)
paire2=btemp(ROI3,ROI4,:)
plot(Time,pair1,'r')
hold on
plot(Time,pair2,'b')

% power spectral for the time serie
PS1= abs(fft(pair1);
PS2= abs(fft(pair2);
plot(PS1,'r',PS2,'b');


% state identification

% exemplar selection
for i=1:N_win
    var=var(Dyn_FC(:,:,i))


% extract features from the ICOV matrix
N_ROI=200
N_comp=7
n=2;
feature=btemp(2,1);
n=n+1:N_ROI
y=[feature, btemp(n,1:n-1)]

% put exemplars from all subjects into a matrix
y_exmp=

% put all windows from all subjects into a matrix
y_total=

% k-mean clustering
[indx,ctrs]=kmeans(y_exemp,N_comp,'Distance','cityblock','Replicates',500)
[indx,ctrs]=kmeans(y_total,N_comp,'Distance','cityblock','start', ctrs)


% calculate percent occurence of each state

for i=1:N_comp
N_occ(i)=sum(indx==i)
perc_occ(i)=N_occ/N_win*N_sub
end

% put the cetroids back as a matix
state1=ctrs(1,:);
state=ones(200,200);
stat

% modularity analysis with Louvain algorithm
for i=1:N_comp
[Ci(i) Q(i)] = modularity_louvain_und_sign(state(i),'sta')
end






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









