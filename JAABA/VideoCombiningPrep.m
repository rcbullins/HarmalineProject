% PURPOSE
%   Preps video organization to be ran on cluster for JAABA training later.
%   Will prep all sessions specified in Directory_Animals Script.
% INPUT
%   Video folder of all raw side and front videos for all trials
% TODO: 
%   Test line 40-45: makes copy of video for you if there is not a
%   copy folder already made. (Previously I hand copied the raw video folders
%   as an extra back up - but the code should now do it for you if a copy does
%   not exist)
% HISTORY
%   1/12/2022: Reagan Bullins
%% Set and add paths
BASEPATH = 'D:\rbullins\'; % Drive and user where raw data and code is located
RAW_DATA = [BASEPATH 'Data\'];% Data location with original files (movies and APT trk files)

JAABAST = [BASEPATH 'Code\JAABAST']; % 
REAGAN_CODE = 'C:\Users\bullinsr\OneDrive - Univeristy of North Carolina at Chapel Hill\Hantman_Lab\Harmaline_Project\Code';
addpath(genpath(REAGAN_CODE));
addpath(genpath(JAABAST));

cd('C:\Users\bullinsr\OneDrive - University of North Carolina at Chapel Hill\Hantman_Lab\Harmaline_Project\Code');
Directory_Animals;
%% Experimental conditions
exper_conditions = {'control';'harm'};
%% For each subject, each condition, each experiment, find the movement trajectory
% Loop through subjects
for isub = 1:length(animals)
    % Identify subject name
    SUB = animals{isub};
    % Loop through each experimental condition
    for iexper = 1:length(exper_conditions)
        % Identify experimental condition
        EXPER_COND = exper_conditions{iexper};
        % Find all sessions in this experimental condition for this subject
        ExperSessions = eval(sprintf('%s_%sBehaviorVideos',SUB,EXPER_COND));
        % Loop through each session
        for isession = 1:length(ExperSessions)
            % Identify experimental session (date of experiment)
            EXPER_SESSION = ExperSessions{isession};
            % If no session exists, skip to next one
            if isempty(EXPER_SESSION)
                continue;
            end
            % Identify path with trk files (estimated position from
            % tracker)
            VIDEO_DIR = [RAW_DATA SUB 'necab1_Chr2/' EXPER_SESSION 'movies/'];
            VIDEO_DIR_COPY = [RAW_DATA SUB 'necab1_Chr2/' EXPER_SESSION 'movies - Copy/'];
            
            if ~exist(VIDEO_DIR_COPY) %have not tested this
            % if video copy does not exist, make a copy
              mkdir(VIDEO_DIR_COPY);
            % Make the new movie and trk file name and directories
              copyfile(VIDEO_DIR, VIDEO_DIR_COPY);
            end
            setUpDir_Part1(VIDEO_DIR_COPY,'frontside',true); 
        end % sessions
    end % conditions
end % animals
%%
%StartJAABA;