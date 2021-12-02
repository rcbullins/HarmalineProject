function [i_cross, data_save] = get_event_ind_v2( fname,nchan,chan2scan,thresh,invert )
% Get indices of digital pulses from binary analog data file. If invert = 1, it
% looks for negative crossings.
% [ i_cross ] = get_event_ind( fname,nchan,chan2scan,thresh,invert )


i_step = 50000*60*5;
%scale = 2.5./(200*(2^15));
f = fopen(fname);
%olap = 50000;
%thresh = .001;
if isempty(thresh)
    thresh = .001*ones(1,length(chan2scan));
end

% Iterate over blocks.
fseek(f,0,'eof');
sz = ftell(f)/(2*nchan);
frewind(f);
n_block = ceil(sz/i_step);
i_cross = cell(1,length(chan2scan));
data_save = [];
for i = 1:n_block
    
    display(['Processing block ' num2str(i) '/' num2str(n_block) '... ']);
    clear data_tmp data_load
    
    % Load data.
    i_offset = (i-1)*i_step;
    for j = 1:length(chan2scan)
        
        data_tmp(j,:) = load_bin_single_channel('int16',fname,nchan,chan2scan(j),i_offset,i_step);
        
        if invert == 1
            data_tmp(j,:) = -data_tmp(j,:);
        end
    end
    data_save = [data_save data_tmp];
    % Find pulses.
    for j = 1:length(chan2scan)
        %{
        [amp_tmp i_cross_tmp] = findpeaks(diff(data_tmp(j,:)),'minpeakheight',thresh(j));
        if ~isempty(i_cross_tmp)
            i_cross{j} = [i_cross{j} i_cross_tmp+(i-1)*i_step];
            amp{j} = amp_tmp;
        end
        %}
        
        i_cross_tmp = find(data_tmp(j,1:(end-1))<thresh(j) & data_tmp(j,2:end)>=thresh(j));
        i_cross{j} = [i_cross{j} i_cross_tmp+(i-1)*i_step];

    end
    
    
end

for i = 1:length(i_cross)
    i_cross{i} = unique(i_cross{i});
end

end