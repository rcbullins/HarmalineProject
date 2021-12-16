function [lfp] = getLFPfromBin(RAW_EPHYS_FILE,BASEPATH, rec_length, SUB, EXPER_SESSION)
% PURPOSE
%   Load lfp file from .bin file. Loads in chunks of one minute of LFP
%   sampled at 30000Hz, and then downsamples to 1250Hz and saves the lfp
%   segment in a mat file. Once all segments are loaded, they are
%   concatenated and saved as a single lfp file.
%
%   Direcotry setup:
%   Stores lfp files in this directory...
%   Basepath > Data_Analyzed > Subject (example: M340) > Neural >
%   SUB_EXPER_SESSION.lfp
%
% INPUTS
%   RAW_EPHYS_FILE   direcotry of .bin file
%   BASEPATH         basepath of computer where data analyzed mat files are saved
%   rec_length       recording length in seconds
%   SUB              mouse number (to locate folder)
%   EXPER_SESSION    date of recording session (to locate folder)
% OUTPUT
%   lfp  struct lfp (channels X time)
% DEPENDENCIES
%   Reagan Code and Britton code
% HISTORY
%    12.15.2021
%%
[EPHYS_PATH, ephys_filename] = fileparts(RAW_EPHYS_FILE);
ephys_filename = [ephys_filename '.bin'];
% Parse the corresponding metafile
meta = ReadMeta(ephys_filename, EPHYS_PATH);
% Get all lfp data
nSamp = floor(rec_length * SampRate(meta));
nSamp_chunk = [1:(60*SampRate(meta)):nSamp];
% Make empty matrix to concat all data
lfp_all = [];
fraction = 0;
fractionVec = [0:1/size(nSamp_chunk,2):1];
f = waitbar(fractionVec(1),['Loading LFP: ' SUB '-' EXPER_SESSION]);
for ichunk = 1:size(nSamp_chunk,2)-1
    lfp_downsampled = [];
    lfp_chunk_file = [BASEPATH 'Data_Analyzed/' SUB '/Neural/tmp_mat_files/' SUB '_' EXPER_SESSION '_lfpChunk' num2str(ichunk) '.mat'];
    if exist(lfp_chunk_file) == 0
        %get waitbar while loading lfp
        waitbar(fractionVec(ichunk),f,['Loading LFP: ' SUB '-' EXPER_SESSION]);
        %chunk data into start and stop sample points
        nSamp_start = nSamp_chunk(1,ichunk);
        % if on last chunk, go to end
        if ichunk == size(nSamp_chunk,2)-1
            nSamp_end = nSamp_chunk(1,ichunk+1);
            numSamp = 60*SampRate(meta);
        else %if not the end
            nSamp_end = nSamp_chunk(1,ichunk+1)-1;
            numSamp = 60*SampRate(meta);
        end
        % dataArray numChan X samples
        lfp_chunk = ReadBin_RB(nSamp_start, numSamp, meta, ephys_filename,EPHYS_PATH);
        chanList = [1:size(lfp_chunk,1)-1];
        if strcmp(meta.typeThis, 'imec')
            lfp_chunk = GainCorrectIM(lfp_chunk, chanList, meta);
        else
            lfp_chunk = GainCorrectNI(lfp_chunk, chanList, meta);
        end
        % downsample lfp from 30000 to 1250
%         ds_factor = length(lfp_chunk)/30000; % lfp once per second
%         ds_factor = ds_factor*1250; % lfp 1250 Hz
        ds_factor = 30000/1250;

        for ichan = 1:length(chanList)
        lfp_downsampled(ichan,:) = lfp_chunk(ichan,1:ds_factor:end);
        end
        save([BASEPATH 'Data_Analyzed/' SUB '/Neural/tmp_mat_files/' SUB '_' EXPER_SESSION '_lfpChunk' num2str(ichunk) '.mat'],'lfp_downsampled');
    else
        load(lfp_chunk_file);
    end
    lfp_all = horzcat(lfp_all ,lfp_downsampled);
end
close(f);
%plot(lfp(ch,:));
lfp = lfp_all;
% delete the chunk files once all saved
for ichunkDel = 1:size(nSamp_chunk,2)-1
    delete([BASEPATH 'Data_Analyzed/' SUB '/Neural/tmp_mat_files/' SUB '_' EXPER_SESSION '_lfpChunk' num2str(ichunk) '.mat']);
end
end
