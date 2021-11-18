function [] = make_jsons_front_side(my_dir)
    
    if ~exist([my_dir '/tracked'],'dir')
        mkdir([my_dir '/tracked'])
    end
    
    % Get list of front and side movies.
    foo = dir([my_dir '/*front*.avi']);
    for i = 1:length(foo)
        fname_vid_front{i} = [foo(i).folder '/' foo(i).name];
        fname_out_front{i} = [my_dir '/tracked/tracked_front_' num2str(i) '.trk'];
    end
    foo = dir([my_dir '/*side*.avi']);
    for i = 1:length(foo)
        fname_vid_side{i} = [foo(i).folder '/' foo(i).name];
        fname_out_side{i} = [my_dir '/tracked/tracked_side_' num2str(i) '.trk'];
    end

    % Encode as JSON files.
    clear data_tmp
    for i = 1:length(fname_vid_front)
        data_tmp(i).movie_files.view_1 = fname_vid_front{i};
        data_tmp(i).output_files.view_1 = fname_out_front{i};
    end
    foo = jsonencode(containers.Map({'toTrack'},{data_tmp}));

    fid = fopen([my_dir '/tracked/info_front.json'], 'w');
    fwrite(fid, foo, 'char');
    fclose(fid);
    
    clear data_tmp
    for i = 1:length(fname_vid_side)
        data_tmp(i).movie_files.view_1 = fname_vid_side{i};
        data_tmp(i).output_files.view_1 = fname_out_side{i};
    end
    foo = jsonencode(containers.Map({'toTrack'},{data_tmp}));

    fid = fopen([my_dir '/tracked/info_side.json'], 'w');
    fwrite(fid, foo, 'char');
    fclose(fid);

end

