function [ out_ts_arr ] = autocorr_resample_ycg( ts, N )
%autocorr_resample generate surrogate data with the same temporal
%   autocorrelation as the input time series

    out_ts_arr = nan;

    % make sure that we are dealing with a single vector and conform
    % it to a column vector
    if min(size(ts)) == 1
        ts=reshape(ts,1,numel(ts));
    else 
        disp(sprintf('Expecting a single time series, but received %d', ...
            min(size(ts))))
    end
    
    out_ts_arr=zeros(length(ts),N);
    
    % perform fourier resampling
    F=fft(ts);
    mag=abs(F);
    phase=angle(F);
    
    
    PhaseRandom_Half = phase(2:length(ts)/2);
    PhaseRandom_Half = PhaseRandom_Half(randperm(length(ts)/2-1));
    %Alternatively, could be: PhaseRandom_Half = rand(1,length(ts)/2-1))*2*pi - pi;
    
    PhaseRandom = [phase(1),PhaseRandom_Half,phase(length(ts)/2+1),-1*flipdim(PhaseRandom_Half,2)];
    
    PhaseRandom = [0,PhaseRandom_Half,0,-1*flipdim(PhaseRandom_Half,2)]; %(0 phase at DC and nyquist freqeuncy)
    
    for i=1:N
        %out_ts_arr(:,i)=mag.*exp(j*phase(randperm(length(ts))));
        out_ts_arr(:,i)=mag.*exp(j*PhaseRandom);
    end
        
    %out_ts_arr=abs(ifft(out_ts_arr));
    out_ts_arr=(ifft(out_ts_arr));
        
end

