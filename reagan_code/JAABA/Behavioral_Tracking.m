%Behavioral_Tracking_Script
% Run this script to prep JAABA
% INPUTS
% OUTPUTS
% DEPENDENCIES 
%   - need .trk files output from APT
%   - need corresponding movie .avi files
%   - need tracker .lbl file from APT tracker
% HISTORY: Reagan Bullins 12/6/2021
%% Add paths
BASEPATH = 'D:/rbullins/'; % Computer at lab only
RAW_DATA = [BASEPATH 'Data/'];
JAABA_CODE = [BASEPATH 'Code/JAABA/perframe'];
JAABA_DATA =  [BASEPATH 'JAABA_behaviors/'];
addpath(genpath(JAABA_DATA));
addpath(JAABA_CODE);

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
            % Identify path with movie files
            MOVIES = [RAW_DATA SUB 'necab1_Chr2/' EXPER_SESSION 'movies/'];
            % Identify path with trk files (estimated position from
            % tracker)
            TRK_FRONT = [RAW_DATA SUB 'necab1_Chr2/' EXPER_SESSION 'movies/trk/front/'];
            TRK_SIDE =  [RAW_DATA SUB 'necab1_Chr2/' EXPER_SESSION 'movies/trk/side/'];
            foldContents = dir(TRK_SIDE);
            numTrials = 0;
            for icont = 1:numel(foldContents)
                foldEmpt = strfind(foldContents(icont).name,'.trk');
                if isempty(foldEmpt)
                    num2add = 0;
                else
                    num2add = 1;
                end
                numTrials = numTrials + num2add;
            end
            % for each trial copy video and trk to new folder for bothe
            % side and front
            for itrial = 1:numTrials
                if itrial < 10
                    thisTrialNum = ['00' num2str(itrial)];
                elseif itrial >=10 && itrial < 100
                    thisTrialNum = ['0' num2str(itrial)];
                elseif itrial >= 100
                    thisTrialNum = num2str(itrial);
                end
                thisTrial_frontMovie = [MOVIES SUB '_' EXPER_SESSION '_front_v' thisTrialNum '.avi'];
                thisTrial_sideMovie = [MOVIES SUB '_' EXPER_SESSION '_side_v' thisTrialNum '.avi'];
                thisTrial_frontTrk = [TRK_FRONT SUB '_' EXPER_SESSION '_front_v' thisTrialNum '.trk'];
                thisTrial_sideTrk = [TRK_SIDE SUB '_' EXPER_SESSION '_side_v' thisTrialNum '.trk'];
                front_fold = [JAABA_DATA SUB '_' EXPER_SESSION '_front_v' thisTrialNum '/'];
                mkdir(front_fold);
                side_fold = [JAABA_DATA SUB '_' EXPER_SESSION '_side_v' thisTrialNum '/'];
                mkdir(side_fold);
                newFront_movie = [JAABA_DATA SUB '_' EXPER_SESSION '_front_v' thisTrialNum '/front/movie.avi'];
                newSide_movie = [JAABA_DATA SUB '_' EXPER_SESSION '_side_v' thisTrialNum '/side/movie.avi'];
                newFront_trk = [JAABA_DATA SUB '_' EXPER_SESSION '_front_v' thisTrialNum '/front/apt_tracks.trk'];
                newSide_trk = [JAABA_DATA SUB '_' EXPER_SESSION '_side_v' thisTrialNum '/side/apt_tracks.trk'];
                copyfile(thisTrial_frontMovie, newFront_movie);
                copyfile(thisTrial_frontTrk, newFront_trk);
                copyfile(thisTrial_sideMovie, newSide_movie);
                copyfile(thisTrial_sideTrk, newSide_trk);
            end

        end %experimental session
    end %experimental condition
end %sub
    
% for each session
% for each trial
% make a new folder with movie copy and trk file copy (with basic file names) if folder does not
% exist

% run JAABA over all folders
StartJAABA;
% Select new for a new behavior
% Select model animal

