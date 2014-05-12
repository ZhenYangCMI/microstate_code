clear
clc

% N: window lenth in TRs
N = 22;
% T: time points
T = 0:N-1;
dr=60;
% sampling rate
fs=0.5;


% creat the rectangular window
RW = ones(1,N);
figure(1); subplot(1,2,1); plot(T,RW)
xlim([0 length(RW)-1])

% Calculate gaussian 
% Gaussian Mean
Gm=(N-1)/2;
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
HanW = hann(N+1)
figure(4); subplot(1,2,1); plot(HanW)
xlim([0 length(HanW)-1])

H4= abs(fftshift(fft(HanW)));
H4= H4/max(H4);
H4= 20*log10(H4);
H4= max(-dr,H4);
h4=-fs/2:fs/(length(H4)):fs/2-fs/(length(H4));
figure(4); subplot(1,2,2);plot(h4,H4)
xlim([-0.25,0.25])

HammingW = hamming(N+1)
figure(5); subplot(1,2,1); plot(HammingW)
xlim([0 length(HammingW)-1])

H5= abs(fftshift(fft(HammingW)));
H5= H5/max(H5);
H5= 20*log10(H5);
H5= max(-dr,H5);
h5=-fs/2:fs/(length(H5)):fs/2-fs/(length(H5));
figure(5); subplot(1,2,2);plot(h5,H5)
xlim([-0.25,0.25])



% cancatenate the time course
TC_seed = load ('Seed_time_series/FunRaw_ROISignals.mat')
TC_atlas= load('Results/FunImgC_ROISignals/ROISignals_0021018.mat')
TC = [TC_seed.theROITimeCourses,TC_atlas.ROISignals];

%stationary FC estimation
cov_matrix = cov(TC)



% invers Fisher transformation


% dynamic FC esimtation and state identification

% segment time course into windows

asize = size(TC)
for i=1:(asize(1,1)-N+1)
    bsize = size(TC(i:N-1+i,:))
    for j= 1:bsize(1,1)
        for k = 1:bsize(1,2)
            w(i,j,k)=TC(i+j-1,k)
        end
    end          
    Cov(i)=cov(w(i,:,:))
end

% create regualrized invese matrix
ts1 = TC
lambda=5; % arbitrary choice of regularisation!
oc=cov(ts1); % raw covariance
ic=-L1precisionBCD(oc/mean(diag(oc)),lambda/1000); % get regularised negative inverse covariance
r=(ic ./ repmat(sqrt(abs(diag(ic))),1,Nnodes)) ./ repmat(sqrt(abs(diag(ic)))',Nnodes,1); % use diagonal to get normalised coefficients
r=r+eye(Nnodes); % remove diagonal 


AM = L1precisionBCD(Cov,1);
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















