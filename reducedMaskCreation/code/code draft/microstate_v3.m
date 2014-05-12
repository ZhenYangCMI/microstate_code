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

% normalize the time series of each voxel
[AllVolume, VoxelSize, ImgFileList, Header, nVolumn] =rest_to4d('swCovRegressed_4DVolume');

for i=1:61
    for j=1:73
        for k=1:61
            for u=1:450
avg=mean(AllVolume(i,j,k,:));
SD=std(AllVolume(i,j,k,:));
z=(AllVolume(i,j,k,u)-avg)/SD;
z(i,j,k,u)=z
            end
        end
    end
end

% 

% cancatenate the normalized time course
TC_seed = load ('data/Seed_time_series/FunRaw_ROISignals.mat')
TC_atlas= load('Results/FunImgC_ROISignals/ROISignals_0021018.mat')
TC = [TC_seed.theROITimeCourses,TC_atlas.ROISignals];


%stationary FC estimation

% if calc_cor=1, generate correlation matrix; if calc_cor=0,calc covariance matrix 
calc_cor = 1

if calc_cor == 1
    stat_FC=corrcoef(TC);
else
    stat_FC=cov(TC);
end
    
% Fisher transformation
z=0.5*log((1+stat_FC)./(1-stat_FC));

% Plot the stationary FC matrix
imagesc(z)
ylabel('seeds and brain ROIs')
xlabel('seeds and brain ROIs')
title('Stationary FC')

% rearrange the regions by lobes
z=z([order of rows],[order of columns]);


% avrage across all subjects



% dynamic FC esimtation 

% Creat matrice containing time series for each seed and the ROIs

for seed=1:4
    TC=[TC(:,seed), TC(:,5:end)]
% segment time course into windows

asize = size(TC);
WC=WC(1,11:32);
N_ROI = 204;
N_win=asize(1,1)-W_width+1;
for i = 1:N_win
    % apply the sliding window to the time series
    for j=2:204
    TC_win_1=WC*TC(i:W_width-1+i,1);
    TC_win=WC*TC(i:W_width-1+i,j);
    TC_win=[TC_win_1,TC_win]
    end
    
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















