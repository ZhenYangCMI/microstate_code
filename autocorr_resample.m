function [ out_ts_arr ] = autocorr_resample( ts, N )
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
        
    for i=1:N
        out_ts_arr(:,i)=mag.*exp(j*phase(randperm(length(ts))));
    end
        
    out_ts_arr=abs(ifft(out_ts_arr));
        
end

