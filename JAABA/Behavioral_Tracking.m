%Behavioral_Tracking_Script
% PURPOSE
%   Run this script to prep data directories for JAABA.
%   This will find all experimental trk files and make copies of them in a
%   directory format that is compatable with JAABA. JAABA requires for
%   each trk file to be in a seperate folder. In this folder (each
%   experimental trial gets its own folder) is the trk file and movie file.
%
%   This script also renames these files to be a general
%   name as necessary for JAABA. After running this script, you may run
%   StartJAABA to load your data in the GUI.
% INPUTS
%   - Change directories in addpaths section
%   - may need to change trk file directories later in script
%       - currently assumes trk file is in ....
%           Sub_ncab1_Chr2/Session_Datemovies/trk/side/
% OUTPUTS
%   - Makes copies of raw trk and movie files, renames to general names
%       apt_tracks.trk and movie.avi, and places them each in their own folder
%       named the original filename (experimental session, front/side, and trial number)
% DEPENDENCIES
%   - .trk files output from APT
%   - corresponding movie .avi files
%   - tracker .lbl file from APT tracker
%   - JAABA downloaded for MATLAB (to run after this script)
%       Download Instructions
%           http://jaaba.sourceforge.net/Installation.html#ComputerRequirements
%       Will direct to clone JAABA from here
%           https://github.com/kristinbranson/JAABA
%       Documentation
%           http://jaaba.sourceforge.net/
% HISTORY: Reagan Bullins 12/6/2021
%% Set and add paths
BASEPATH = 'D:/rbullins/'; % Drive and user where raw data and code is located
RAW_DATA = [BASEPATH 'Data/'];% Data location with original files (movies and APT trk files)
JAABA_DATA =  [BASEPATH 'JAABA_behaviors/'];% NEW path to save copies of files to

JAABA_CODE = [BASEPATH 'Code/JAABA/perframe']; % path of JAABA code cloned from github
REAGAN_CODE = 'C:/Users/bullinsr/OneDrive - Univeristy of North Carolina at Chapel Hill/Hantman_Lab/Harmaline_Project/Code/JAABA';
REAGAN_DIRECTORY = 'C:\Users\bullinsr\OneDrive - University of North Carolina at Chapel Hill\Hantman_Lab\Harmaline_Project\Code';

addpath(REAGAN_DIRECTORY);
addpath(genpath(JAABA_DATA));
addpath(JAABA_CODE);
addpath(genpath(REAGAN_CODE));

%% Load in experimental information of all experimental sessions
Directory_Animals;
%% Experimental conditions
exper_conditions = {'control';'harm'};
%%
% Set directory
% for each animal
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
            % Identify path with movie files of trials
            MOVIES = [RAW_DATA SUB 'necab1_Chr2/' EXPER_SESSION 'movies/'];
            % Identify path with trk files (estimated position from
            % tracker)
            TRK_FRONT = [RAW_DATA SUB 'necab1_Chr2/' EXPER_SESSION 'movies/trk/front/'];
            TRK_SIDE =  [RAW_DATA SUB 'necab1_Chr2/' EXPER_SESSION 'movies/trk/side/'];
            foldContents = dir(TRK_SIDE);
            numTrials = 0;
            % Find how many trials there are with videos and trk files
            for icont = 1:numel(foldContents)
                foldEmpt = strfind(foldContents(icont).name,'.trk');
                if isempty(foldEmpt)
                    num2add = 0;
                else
                    num2add = 1;
                end
                numTrials = numTrials + num2add;
            end
            % for each trial copy video and trk to new folder for both
            % side and front
            for itrial = 1:numTrials
                % get trial number in string form
                if itrial < 10
                    thisTrialNum = ['00' num2str(itrial)];
                elseif itrial >=10 && itrial < 100
                    thisTrialNum = ['0' num2str(itrial)];
                elseif itrial >= 100
                    thisTrialNum = num2str(itrial);
                end
                % Identify movie and trk file locations for this trial
                % session
                thisTrial_frontMovie = [MOVIES SUB '_' EXPER_SESSION '_front_v' thisTrialNum '.avi'];
                thisTrial_sideMovie = [MOVIES SUB '_' EXPER_SESSION '_side_v' thisTrialNum '.avi'];
                thisTrial_frontTrk = [TRK_FRONT SUB '_' EXPER_SESSION '_front_v' thisTrialNum '.trk'];
                thisTrial_sideTrk = [TRK_SIDE SUB '_' EXPER_SESSION '_side_v' thisTrialNum '.trk'];
                % Identify and make folders of where to store the copies
                exper_fold = [JAABA_DATA SUB '_' EXPER_SESSION '_v' thisTrialNum '/'];
                if exist(exper_fold) ==0
                    mkdir(exper_fold);

                    % Make the new movie and trk file name and directories
                    newFront_movie = [JAABA_DATA SUB '_' EXPER_SESSION '_v' thisTrialNum '/movie1.avi'];
                    newSide_movie = [JAABA_DATA SUB '_' EXPER_SESSION '_v' thisTrialNum '/movie2.avi'];
                    newFront_trk = [JAABA_DATA SUB '_' EXPER_SESSION '_v' thisTrialNum '/apt_tracks1.trk'];
                    newSide_trk = [JAABA_DATA SUB '_' EXPER_SESSION '_v' thisTrialNum '/apt_tracks2.trk'];
                    % Copy the files to the new directories with the new name
                    copyfile(thisTrial_frontMovie, newFront_movie);
                    copyfile(thisTrial_frontTrk, newFront_trk);
                    copyfile(thisTrial_sideMovie, newSide_movie);
                    copyfile(thisTrial_sideTrk, newSide_trk);
                end
            end %trials
        end %experimental session
    end %experimental condition
end %sub

%% To Run JAABA
%StartJAABA;


