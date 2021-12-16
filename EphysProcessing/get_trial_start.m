function [tstamp]=get_trial_start(ind_ctx,index)
% HISTORY
%   From Britton 11.23.2021
    X=ind_ctx.ind{index};  
    tstamp=[];
    for i=1:length(X)
        if i==1
            tstamp=[tstamp X(i)];
        else

            X0=X(i-1);
            X1=X(i);
            diff=X1-X0;
            if diff>1000
                tstamp=[tstamp X1];
            end
        end
    end

end



