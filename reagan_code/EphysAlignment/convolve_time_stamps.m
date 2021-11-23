function [rates,t_rates]=convolve_time_stamps(st_ctx,window,sess_len,causal_flag)
% HISTORY
%   From Britton
    nNeurons=size(st_ctx,2);
    bin=1;
    clear rates
    for n=1:nNeurons
       tmp_spks=st_ctx{n};
       train = histc(tmp_spks,0:bin:sess_len);
       kernel = 1000*(1/(window*sqrt(2*pi)))*exp((-1/(2*(window^2)))*((-(window*5):bin:(window*5)).^2));
       if causal_flag
%            kernel=kernel(round(length(kernel))/2:end);
           kernel=kernel(1:round(length(kernel))/2);
       end
       rates(n,:)=conv(train,kernel,'same');
    end
    t_rates=0:1:sess_len; %in ms

end