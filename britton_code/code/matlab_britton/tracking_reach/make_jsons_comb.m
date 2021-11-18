function [] = make_jsons_comb(my_dir,frame_range)

    % Set parameters ...
    %my_dir_UNIX = ['/groups/hantman/' my_dir(10:end)];
    my_dir_UNIX = my_dir;
    if nargin == 1
        frame_range = [1 2501];
    end
    mkdir([my_dir '/tracked_fingers']);

    % Get list of combined movies.
    %{
    foo = dir([my_dir '/*/*/movie_comb.avi']);
    if isempty(foo)
        foo = dir([my_dir '/*/movie_comb.avi']);
    end
    %}
    foo = dir([my_dir '/*/movie_comb.avi']);
    
    % SORT BY TRIAL NUMBER!
    trial = [];
    for i = 1:length(foo)
        bar = strfind(foo(i).folder,'_v');
        trial(i) = str2num(foo(i).folder(bar+(2:4)));
    end
    [trial order] = sort(trial);
    foo = foo(order);
    
    for i = 1:length(foo)
        fname_vid_comb{i} = [foo(i).folder '/' foo(i).name];
        frames_comb{i} = {i 1 frame_range};
    end

    % Encode as JSON file.
    foo = jsonencode(containers.Map({'movieFiles' 'toTrack'},...
    {fname_vid_comb frames_comb}));
    fid = fopen([my_dir '/tracked_fingers/info_comb.json'], 'w');
    fwrite(fid, foo, 'char');
    fclose(fid);


    % Save info.
    jsons = {[my_dir_UNIX '/tracked_fingers/info_comb.json']};
    outfiles = {[my_dir_UNIX '/tracked_fingers/tracked_comb']};
    save([my_dir '/tracked_fingers/job_info.mat'],'jsons','outfiles','trial')

end


