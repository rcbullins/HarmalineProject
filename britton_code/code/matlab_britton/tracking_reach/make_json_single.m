function [] = make_json_single(vid_name,frame_range,obj2track,out_name)
% function [] = make_json_single(vid_name,frame_range,obj2track,out_name)

    % Set parameters ...
    foo = strfind(vid_name,'/');
    my_dir = vid_name(1:(foo(end)-1));
    my_dir_UNIX = ['/groups/hantman/' my_dir(10:end)];
    vid_name_UNIX = [my_dir_UNIX vid_name(foo(end):end)];
    foo = strfind(out_name,'/');
    out_name_short = out_name((foo(end)+1):end);
    
    for i = 1:length(obj2track)
        ii = obj2track(i);
        my_vid{i} = vid_name_UNIX;
        my_frames{i} = {1 ii frame_range};
    end
    
    % Encode as JSON files.
    foo = jsonencode(containers.Map({'movieFiles' 'toTrack'},...
    {my_vid my_frames}));
    fid = fopen([my_dir '/info_' out_name_short '.json'], 'w');
    fwrite(fid, foo, 'char');
    fclose(fid);

    % Save info.
    jsons = {[my_dir_UNIX '/info_' out_name_short '.json']};
    outfiles = {[my_dir_UNIX '/' out_name_short]};
    save([my_dir '/job_info.mat'],'jsons','outfiles')

end

