function [] = make_jsons_front_side_ALL(dirs_all,dir_tmp)
    
    fname_vid_front = [];
    fname_out_front = [];
    fname_vid_side = [];
    fname_out_side = [];
    % Loop over datasets and get all the video and output locations.
    for dd = 1:length(dirs_all)
        if ~exist([dirs_all{dd} '/tracked'],'dir')
            mkdir([dirs_all{dd} '/tracked'])
        end

        % Get list of front and side movies.
        i_start = length(fname_vid_front);
        foo = dir([dirs_all{dd} '/*front*.avi']);
        for i = 1:length(foo)
            fname_vid_front{i_start+i} = [foo(i).folder '/' foo(i).name];
            fname_out_front{i_start+i} = [dirs_all{dd} '/tracked/tracked_front_' num2str(i) '.trk'];
        end
        
        i_start = length(fname_vid_side);
        foo = dir([dirs_all{dd} '/*side*.avi']);
        for i = 1:length(foo)
            fname_vid_side{i_start+i} = [foo(i).folder '/' foo(i).name];
            fname_out_side{i_start+i} = [dirs_all{dd} '/tracked/tracked_side_' num2str(i) '.trk'];
        end
    end

    % Encode as JSON files.
    clear data_tmp
    for i = 1:length(fname_vid_front)
        data_tmp(i).movie_files.view_1 = fname_vid_front{i};
        data_tmp(i).output_files.view_1 = fname_out_front{i};
    end
    foo = jsonencode(containers.Map({'toTrack'},{data_tmp}));

    fid = fopen([dir_tmp '/info_front.json'], 'w');
    fwrite(fid, foo, 'char');
    fclose(fid);
    display(['Front view *.json written, ' num2str(length(fname_vid_front)) ' videos.']);
    
    clear data_tmp
    for i = 1:length(fname_vid_side)
        data_tmp(i).movie_files.view_1 = fname_vid_side{i};
        data_tmp(i).output_files.view_1 = fname_out_side{i};
    end
    foo = jsonencode(containers.Map({'toTrack'},{data_tmp}));

    fid = fopen([dir_tmp '/info_side.json'], 'w');
    fwrite(fid, foo, 'char');
    fclose(fid);
    display(['Side view *.json written, ' num2str(length(fname_vid_side)) ' videos.']);
    
    if length(fname_vid_front) ~= length(fname_vid_side)
        display('WARNING: mismatch between numbers of front and side view videos');
    end
    
end

