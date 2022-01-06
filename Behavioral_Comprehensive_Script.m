% Behavior script Harmaline Project
% PURPOSE
%   Run all necessary behavior analysis for harmaline dataset. This
%   script will create 2D and 3D trajectories of reach to grab movements,
%   calculate accuracy, look at endpoints, compare control to harmaline experimental
%   sessions.
% DEPENDENCIES
%   Britton Code
%       Animal Tracking Calibration Files
%       Original get_traj_3D code
%   Trk files
%       Output from APT tracking (estimated x,y,z position of body parts)
%           - One trk file per trial (reach to grab trial)
%   JAABA output
%       JAABA is a more accurate way to define parts of the reach to grab
%       movement. See Behavioral_Tracking script to run JAABA.
%   Excel experimental session information sheets
%       One for each subject
%       Created at the time of experiment by Jay
%       Each session: clear background = baseline, colored background = stim
%       General session sequence: baseline > stim VA ChrRhp > washout
%   Directory_Animals.m
%       Script outlining experimental info: sessions and subjects to run
%       analysis over. Consult excel sheets for identifying trials as
%       baseline, stimulation, or washout. (already done, but to add
%       sessions edit this script.) NOTE: adding more than one session to
%       each condition (control vs harmaline) might require some editing of
%       code to account for that. :)
%   Reagan_code for harmaline project
%       required functions to run script
%       getStartFrame (of movement)
%       getEndFrame (of movement)
%       plot1DTrajectories
%       plot3DTrajectories
%       plotOverlay1D
%       plotOverlay3D
%       Behavioral_Comparison_Conditions (plots 3D trajectories comparing
%       harmaline and control subjects, plots endpoints only, aligns
%       endpoints based on pellet location)
%       Behavioral_Metrics (plots accuracy comparing control vs harmaline
%       and baseline vs stimulation trials)
% OUTPUTS
%       See section below to specify what graphs to output: Samples
%           plotSampleTrajectories (plot 2D and 3D sample trajectories)
%           plot1DOverlayTrajectories (plot all 2D trajectories on same
%               graph)
%           plot3DOverlayTrajectories (plot all 3D trajectories on same
%               graph)
%           plotAccuracy(plot individual accuracy for each condition)
%       Comparison graphs (will always make these)
%           Comparing accuracy on reach movement between harmaline and
%               control subject and on baseline vs stimulation trials
%           Comparing 3D trajectories (from movement initation to pellet
%               grab) of harmaline vs control dataset
%               - Another plotting endpoints only (harmaline red, control
%                    blue) Dots are pellet locations
%               - Another aligning control and hamrlaine sessions using the
%                    location of the pellet. Then plotting endpoints as success
%                    (o) or failures (x)
% HISTORY
%   11/23/2021 Reagan Bullins
%% Clean workspace
clear;
close ALL;
clc;
%% Specify what to plot
useJAABA = 1;
plotSampleTrajectories = 0;
plot1DOverlayTrajectories = 0;
plot3DOverlayTrajectories = 0;
plotAccuracy = 0;
score = 'all'; % Options: 1, 0, 2, -1, 'all'
% Code: (1) one grab and success
%       (0) grab and failure
%       (2) multiple reaches and eventual success
%      (-1) no reach attempts
%   ('all') all scores where some attempt was made
%% Add code paths
USER = 'bullinsr';
RAWDATA_BASEPATH = 'D:/rbullins/'; % Computer at lab only
BASEPATH = ['C:/Users/' USER '/OneDrive - University of North Carolina at Chapel Hill/Hantman_Lab/Harmaline_Project/'];
CODE_REAGAN = [BASEPATH 'Code/reagan_code/'];
CODE_CALIB = [RAWDATA_BASEPATH 'Code/britton_code/calibration_files/camera_calibration_Jay_7-28-16/Calib_Results_stereo.mat'];
CODE_TRKR = [RAWDATA_BASEPATH 'Code/britton_code/tracker'];
CODE_BRITTON = [RAWDATA_BASEPATH 'Code/britton_code/code/matlab_britton/'];
CODE_BRITTON_PLOT = [RAWDATA_BASEPATH 'Code/britton_code/other_code'];
CODE_PROCESS_EVENTS = [RAWDATA_BASEPATH 'Code/process_events/'];
HARM_CODE = [BASEPATH 'Code/'];
JAABA_OUTPUT = [RAWDATA_BASEPATH 'JAABA_behaviors/'];
addpath(genpath(JAABA_OUTPUT));
addpath(genpath(HARM_CODE));
addpath(genpath(CODE_REAGAN));
addpath(CODE_CALIB);
addpath(CODE_TRKR);
addpath(genpath(CODE_BRITTON));
addpath(genpath(CODE_BRITTON_PLOT));
addpath(genpath(CODE_PROCESS_EVENTS));
% Add paths to data inventory script (edit to specify sessions to run over)
Directory_Animals;
%% Experimental conditions
exper_conditions = {'control';'harm'};
%% Score
scoreLabel = num2str(score);
if strcmp(scoreLabel, '1')
    SCORE = 'idealSuccess';
elseif strcmp(scoreLabel, '0')
    SCORE = 'noSuccess';
elseif strcmp(scoreLabel, '-1')
    SCORE = 'noReach';
elseif strcmp(scoreLabel, '2')
    SCORE = 'eventualSuccess';
elseif strcmp(scoreLabel, 'all')
    score = [1 0 2];
    SCORE = 'allScores';
end
%% For each subject, each condition, each experiment, find the movement trajectory
% Loop through subjects
for isub = 1 %:length(animals)
    % Identify subject name
    SUB = animals{isub};
    % Loop through each experimental condition
    for iexper = 1%:length(exper_conditions)
        % Identify experimental condition
        EXPER_COND = exper_conditions{iexper};
        % Find all sessions in this experimental condition for this subject
        ExperSessions = eval(sprintf('%s_%sBehaviorVideos',SUB,EXPER_COND));
        % Add path for where to save figures
        BEHAV_FIG = [BASEPATH 'Figures/' SUB '/Behavior/'];
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
            TRK = [RAWDATA_BASEPATH 'Data/' SUB 'necab1_Chr2/' EXPER_SESSION 'movies/trk'];
            addpath(TRK);
            %% Load position data
            % For digit 1
            trajFile = ([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_trajectories.mat']);
            if exist(trajFile,'file')==2
                load(trajFile);
            else
                nframes=2500; %variable! %number of frames to analyze. Going to depend on how many frames in movie
                digit2trk = 1;
                [traj, ~] = get_traj_3D_v2(TRK,CODE_CALIB,nframes,digit2trk);
                save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_trajectories.mat'],'traj');
            end
            % For digit 2
            trajFileTwo = ([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_trajectoriesDigitTwo.mat']);
            if exist(trajFileTwo,'file')==2
                load(trajFileTwo);
            else
                digit2trk = 2;
                nframes = 2500;
                [trajDigTwo, ~] = get_traj_3D_v2(TRK,CODE_CALIB,nframes,digit2trk);
                save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_trajectoriesDigitTwo.mat'],'trajDigTwo');
            end
            % traj (trials x coordinate xyz x frames);
            % conf is the confidence value assigned by the classifier
            %% Set trials
            % Define Trials by baseline, stim, washout
            trialIdxs = eval(sprintf('%s_%s_%sTrials',SUB,EXPER_SESSION,EXPER_COND));
            base.trialIdxs = trialIdxs.nbase;
            pert.trialIdxs = trialIdxs.npert;
            wash.trialIdxs = trialIdxs.nwash;
            %% Assign trajectories based on trial type
            % traj(trials x coordinate xyz x frames);
            base.traj = traj(base.trialIdxs,:,:);
            pert.traj = traj(pert.trialIdxs,:,:);
            wash.traj = traj(wash.trialIdxs,:,:);
            %% Set Graph Defaults - I like arial font :)
            SetGraphDefaults;
            %% Get num of success, no success, etc for each baseline,stim, and wash
            % baseline
            numBL.idealSuccess = length(find(trialIdxs.trialScore(base.trialIdxs) == 1));
            numBL.eventualSuccess = length(find(trialIdxs.trialScore(base.trialIdxs) == 2));
            numBL.noSuccess = length(find(trialIdxs.trialScore(base.trialIdxs) == 0));
            numBL.noReach = length(find(trialIdxs.trialScore(base.trialIdxs) == -1));
            numBL.grooming = length(find(trialIdxs.trialScore(base.trialIdxs) == 'g'));
            % stim
            numPert.idealSuccess = length(find(trialIdxs.trialScore(pert.trialIdxs)==1));
            numPert.eventualSuccess = length(find(trialIdxs.trialScore(pert.trialIdxs) == 2));
            numPert.noSuccess = length(find(trialIdxs.trialScore(pert.trialIdxs) == 0));
            numPert.noReach = length(find(trialIdxs.trialScore(pert.trialIdxs) == -1));
            numPert.grooming = length(find(trialIdxs.trialScore(pert.trialIdxs) == 'g'));
            % washout
            numWash.idealSuccess = length(find(trialIdxs.trialScore(wash.trialIdxs)==1));
            numWash.eventualSuccess = length(find(trialIdxs.trialScore(wash.trialIdxs) == 2));
            numWash.noSuccess = length(find(trialIdxs.trialScore(wash.trialIdxs) == 0));
            numWash.noReach = length(find(trialIdxs.trialScore(wash.trialIdxs) == -1));
            numWash.grooming = length(find(trialIdxs.trialScore(wash.trialIdxs) == 'g'));
            %% Find accuracy for baseline and save mat
            % Find total number of baseline trials
            numTotalBL = length(base.trialIdxs);
            % Find number of each trial type
            y_BL = [numBL.idealSuccess; numBL.eventualSuccess ;numBL.noSuccess ;numBL.noReach; numBL.grooming];
            % Divide each trial type by the total number to get percentage
            y_BL_acc = y_BL./numTotalBL;
            % Save the accuracy as a mat file to be called upon later in
            % Behavioral_Metrics
            save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_baselineAccuracy.mat'],'y_BL_acc');
            % Plot Accuracy (will be plotted in conjunction with other
            % experimental conditions later)
            if plotAccuracy == 1
                h = bar(y_BL_acc*100,'FaceColor','flat');
                colorScheme = colorGradient([1,0,.5],[.5,0,1],5);
                h.CData(1,:) = colorScheme(1,:);
                h.CData(2,:) = colorScheme(2,:);
                h.CData(3,:) = colorScheme(3,:);
                h.CData(4,:) = colorScheme(4,:);
                h.CData(5,:) = colorScheme(5,:);
                l = cell(1,5);
                l{1}='Success'; l{2}='Eventual Success'; l{3}='No Success'; l{4}='No Reach'; l{5}='Grooming';
                set(gca,'xticklabel', l);
                ylabel('Accuracy (%)');
                title([SUB ' ' EXPER_SESSION ' ' EXPER_COND ' :BL Accuracy']);
            end
            %% Find stimulation accruacy and save mat
            numTotalPert = length(pert.trialIdxs);
            y_Pert = [numPert.idealSuccess; numPert.eventualSuccess ;numPert.noSuccess ;numPert.noReach; numPert.grooming];
            y_Pert_acc = y_Pert./numTotalPert;
            save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_pertAccuracy.mat'],'y_Pert_acc');
            if plotAccuracy == 1
                figure;
                h = bar(y_Pert_acc*100,'FaceColor','flat');
                colorScheme = colorGradient([1,0,.5],[.5,0,1],5);
                h.CData(1,:) = colorScheme(1,:);
                h.CData(2,:) = colorScheme(2,:);
                h.CData(3,:) = colorScheme(3,:);
                h.CData(4,:) = colorScheme(4,:);
                h.CData(5,:) = colorScheme(5,:);
                l = cell(1,5);
                l{1}='Success'; l{2}='Eventual Success'; l{3}='No Success'; l{4}='No Reach'; l{5}='Grooming';
                set(gca,'xticklabel', l);
                ylabel('Accuracy (%)');
                title([SUB ' ' EXPER_SESSION ' ' EXPER_COND ' : Stim Accuracy']);
                %            savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajector/' SUB '_' EXPER_SESSION '_' EXPER_COND '_zDim_washout.fig']);
            end
            %% Find accuracy of eventual success /only failed reaches and multiple reaches
            % Stimulation
            numTotalIsolatePert = numPert.eventualSuccess + numPert.noSuccess;
            y_Pert_isolate = [numPert.eventualSuccess];
            y_Pert_acc_isolate = y_Pert_isolate./numTotalIsolatePert;
            save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_pertIsolateAccuracy.mat'],'y_Pert_acc_isolate');
            % Baseline
            numTotalIsolateBL = numBL.eventualSuccess + numBL.noSuccess;
            y_BL_isolate = [numBL.eventualSuccess];
            y_BL_acc_isolate = y_BL_isolate./numTotalIsolateBL;
            save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_BLIsolateAccuracy.mat'],'y_BL_acc_isolate');

            %% Plot Trajectory examples (first 16 trials)
            if plotSampleTrajectories == 1
                %% Plot Sample 1D Trajectories in Z
                nframes = 2500;
                % INPUT PARSERS AVAILABLE: coordinate_dim, startFrame, endFrame
                plot1DTrajectories(base.traj,'endFrame',nframes);
                sgtitle([SUB ' ' EXPER_COND ': handpaths baseline trials ZDim']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_zDim_baseline.fig']);
                plot1DTrajectories(pert.traj, 'endFrame',nframes);
                sgtitle([SUB ' ' EXPER_COND ': handpaths stimulation trials ZDim']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_zDim_stim.fig']);
                plot1DTrajectories(wash.traj, 'endFrame',nframes);
                sgtitle([SUB ' ' EXPER_COND ': handpaths washout trials ZDim']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_zDim_washout.fig']);
                %% Plot Sample 1D Trajectories in y
                % INPUT PARSERS AVAILABLE: coordinate_dim, startFrame, endFrame
                plot1DTrajectories(base.traj,'endFrame',nframes,'coordinate_dim',2);
                sgtitle([SUB ' ' EXPER_COND ': handpaths baseline trials YDim']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_yDim_baseline.fig']);
                plot1DTrajectories(pert.traj, 'endFrame',nframes,'coordinate_dim',2);
                sgtitle([SUB ' ' EXPER_COND ': handpaths stimulation trials YDim']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_yDim_stim.fig']);
                plot1DTrajectories(wash.traj, 'endFrame',nframes,'coordinate_dim',2);
                sgtitle([SUB ' ' EXPER_COND ': handpaths washout trials YDim']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_yDim_washout.fig']);

                %% Plot 3D Trajectories
                plot3DTrajectories(base.traj, 'endFrame',nframes);
                sgtitle([SUB ' ' EXPER_COND ': 3D handpaths baseline trials']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/3DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_baseline.fig']);
                plot3DTrajectories(pert.traj, 'endFrame',nframes);
                sgtitle([SUB ' ' EXPER_COND ': 3D handpaths stimulation trials']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/3DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_stim.fig']);
                plot3DTrajectories(wash.traj, 'endFrame',nframes);
                sgtitle([SUB ' ' EXPER_COND ': 3D handpaths washout trials']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/3DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_washout.fig']);
            end
            %% Find start frame of movement (uses smoothing and thresholding of velocity)
            % Only for trials with a score of 'score' variable
            [movStrt.nbase] = getMovStartFrame(score,traj,trialIdxs.trialScore,trialIdxs.nbase);
            [movStrt.npert] = getMovStartFrame(score,traj,trialIdxs.trialScore,trialIdxs.npert);
            [movStrt.nwash] = getMovStartFrame(score,traj,trialIdxs.trialScore,trialIdxs.nwash);
            %% Find end frame of reach movement (uses smoothing and thresholding of velocity)
            % Only for trials with a score of 'score' defined at the start
            % of the script.
            % movIdx is the index of the trial within all trials
            [movEnd.nbase,movIdx.nbase] = getMovEndFrame(score,traj,trialIdxs.trialScore,trialIdxs.nbase, movStrt.nbase);
            [movEnd.npert,movIdx.npert] = getMovEndFrame(score,traj,trialIdxs.trialScore,trialIdxs.npert, movStrt.npert);
            [movEnd.nwash,movIdx.nwash] = getMovEndFrame(score,traj,trialIdxs.trialScore,trialIdxs.nwash, movStrt.nwash);
            %% Using JAABA - Start and End Frames of reach and grab
            if useJAABA
                % Start: Hand Lift and End: Grab
                % Load classifier data (hand lift, grab, etc)
                % find how many trials there are
                TRK_SIDE =  [RAWDATA_BASEPATH 'Data/' SUB 'necab1_Chr2/' EXPER_SESSION 'movies/trk/side/'];
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
                % Initiate cells to store start frames in
                chew.t0           = {};
                atMouth.t0        = {};
                digitsTogether.t0 = {};
                handLift.t0       = {};
                supinate.t0       = {};
                grab.t0           = {};
                % Initiate cells to store stop frames in
                chew.t1           = {};
                atMouth.t1        = {};
                digitsTogether.t1 = {};
                handLift.t1       = {};
                supinate.t1       = {};
                grab.t1           = {};
                % load output for all trials and store in mega matrix
                for itrial = 1:numTrials
                    if itrial < 10
                        thisTrialNum = ['00' num2str(itrial)];
                    elseif itrial >=10 && itrial < 100
                        thisTrialNum = ['0' num2str(itrial)];
                    elseif itrial >= 100
                        thisTrialNum = num2str(itrial);
                    end
                    TRIAL_FOLDER = [JAABA_OUTPUT SUB '_' EXPER_SESSION '_v' thisTrialNum '/'];
                    chew_scores           = load([TRIAL_FOLDER 'scores_Chew.mat'],'allScores');
                    atMouth_scores        = load([TRIAL_FOLDER 'scores_AtMouth.mat'],'allScores');
                    digitsTogether_scores = load([TRIAL_FOLDER 'scores_DigitsTogether.mat'],'allScores');
                    handLift_scores       = load([TRIAL_FOLDER 'scores_LiftHand.mat'],'allScores');
                    supinate_scores       = load([TRIAL_FOLDER 'scores_Supinate.mat'],'allScores');
                    grab_scores           = load([TRIAL_FOLDER 'scores_Grab.mat'],'allScores');
                    % save in mega cell array start times
                    chew.t0{itrial}           = chew_scores.allScores.t0s;
                    atMouth.t0{itrial}        = atMouth_scores.allScores.t0s;
                    digitsTogether.t0{itrial} = digitsTogether_scores.allScores.t0s;
                    handLift.t0{itrial}       = handLift_scores.allScores.t0s;
                    supinate.t0{itrial}       = supinate_scores.allScores.t0s;
                    grab.t0{itrial}           = grab_scores.allScores.t0s;
                    % save in mega cell array stop times
                    chew.t1{itrial}           = chew_scores.allScores.t1s;
                    atMouth.t1{itrial}        = atMouth_scores.allScores.t1s;
                    digitsTogether.t1{itrial} = digitsTogether_scores.allScores.t1s;
                    handLift.t1{itrial}       = handLift_scores.allScores.t1s;
                    supinate.t1{itrial}       = supinate_scores.allScores.t1s;
                    grab.t1{itrial}           = grab_scores.allScores.t1s;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % find move start based on hand lift with indexes for each
                % BASELINE start frame and end frame
                movStart.nbase = zeros(1,length(movIdx.nbase));
                movEnd.nbase = zeros(1,length(movIdx.nbase));
                for itrial = 1:length(movIdx.nbase) %1:numTrials
                    % assume start is first hand lift
                    if isempty(handLift.t0{1,movIdx.nbase(itrial)}{1,1}) == 0
                        movStart.nbase(1,itrial) = handLift.t0{1,movIdx.nbase(itrial)}{1,1}(1);
                        noLift.nbase = 0;
                    else
                        movStart.nbase(1,itrial) = 0;
                        noLift.nbase = 1;
                    end
                    % assume end is first grab
                    if isempty(grab.t0{1,movIdx.nbase(itrial)}{1,1}) == 0
                        movEnd.nbase(1,itrial) = grab.t0{1,movIdx.nbase(itrial)}{1,1}(1);
                        noGrab.nbase = 0;
                    else
                        movEnd.nbase(1,itrial) = 0;
                        noGrab.nbase = 1;
                    end
                    % if there is neither lift or grab, move to next trial
                    if noLift.nbase || noGrab.nbase
                        continue;
                    end
                    % Now lets check to see that grab happens after lift
                    if movStart.nbase(1,itrial) < movEnd.nbase(1,itrial)
                        % this is good, let's check there is not a lift
                        % between this lift and the grab.
                        liftsBeforeGrab = find(handLift.t0{1,movIdx.nbase(itrial)}{1,1} < movEnd.nbase(1,itrial));
                        if liftsBeforeGrab > 1
                            %this means there are more than one lifts
                            %before the grab, we should choose the one
                            %closest to the grab
                            liftClosest2Grab = handLift.t0{1,movIdx.nbase(itrial)}{1,1}(length(liftsBeforeGrab));
                            movStart.nbase(1,itrial) = liftClosest2Grab;
                        end % if it is one, then that is the only lift
                    else % if the grab happens before the lift, this is an early grab. Find the grab after the lift.
                        % find grab after lift and after digits together
                        firstGrabAfterLift = find(grab.t0{1,movIdx.nbase(itrial)}{1,1} > handLift.t0{1,movIdx.nbase(itrial)}{1,1}(1));

                        movEnd.nbase(1,itrial) = grab.t0{1,movIdx.nbase(itrial)}{1,1}(firstGrabAfterLift(1));

                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % PERT start frame and end frame
                movStart.npert = zeros(1,length(movIdx.npert));
                movEnd.npert = zeros(1,length(movIdx.npert));
                for itrial = 1:length(movIdx.npert) %1:numTrials
                    % assume start is first hand lift
                    if isempty(handLift.t0{1,movIdx.npert(itrial)}{1,1}) == 0
                        movStart.npert(1,itrial) = handLift.t0{1,movIdx.npert(itrial)}{1,1}(1);
                        noLift.npert = 0;
                    else
                        movStart.npert(1,itrial) = 0;
                        noLift.npert = 1;
                    end
                    % assume end is first grab
                    if isempty(grab.t0{1,movIdx.npert(itrial)}{1,1}) == 0
                        movEnd.npert(1,itrial) = grab.t0{1,movIdx.npert(itrial)}{1,1}(1);
                        noGrab.npert = 0;
                    else
                        movEnd.npert(1,itrial) = 0;
                        noGrab.npert = 1;
                    end
                    % if there is neither lift or grab, move to next trial
                    if noLift.npert || noGrab.npert
                        continue;
                    end
                    % Now lets check to see that grab happens after lift
                    if movStart.npert(1,itrial) < movEnd.npert(1,itrial)
                        % this is good, let's check there is not a lift
                        % between this lift and the grab.
                        liftsBeforeGrab = find(handLift.t0{1,movIdx.npert(itrial)}{1,1} < movEnd.npert(1,itrial));
                        if liftsBeforeGrab > 1
                            %this means there are more than one lifts
                            %before the grab, we should choose the one
                            %closest to the grab
                            liftClosest2Grab = handLift.t0{1,movIdx.npert(itrial)}{1,1}(length(liftsBeforeGrab));
                            movStart.npert(1,itrial) = liftClosest2Grab;
                        end % if it is one, then that is the only lift
                    else % if the grab happens before the lift, this is an early grab. Find the grab after the lift.
                        % find grab after lift and after digits together
                        firstGrabAfterLift = find(grab.t0{1,movIdx.npert(itrial)}{1,1} > handLift.t0{1,movIdx.npert(itrial)}{1,1}(1));

                        movEnd.npert(1,itrial) = grab.t0{1,movIdx.npert(itrial)}{1,1}(firstGrabAfterLift(1));
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %WASHOUT start frame and end frame
                movStart.nwash = zeros(1,length(movIdx.nwash));
                movEnd.nwash = zeros(1,length(movIdx.nwash));
                for itrial = 1:length(movIdx.nwash) %1:numTrials
                    % assume start is first hand lift
                    if isempty(handLift.t0{1,movIdx.nwash(itrial)}{1,1}) == 0
                        movStart.nwash(1,itrial) = handLift.t0{1,movIdx.nwash(itrial)}{1,1}(1);
                        noLift.nwash = 0;
                    else
                        movStart.nwash(1,itrial) = 0;
                        noLift.nwash = 1;
                    end
                    % assume end is first grab
                    if isempty(grab.t0{1,movIdx.nwash(itrial)}{1,1}) == 0
                        movEnd.nwash(1,itrial) = grab.t0{1,movIdx.nwash(itrial)}{1,1}(1);
                        noGrab.nwash = 0;
                    else
                        movEnd.nwash(1,itrial) = 0;
                        noGrab.nwash = 1;
                    end
                    % if there is neither lift or grab, move to next trial
                    if noLift.nwash || noGrab.nwash
                        continue;
                    end
                    % Now lets check to see that grab happens after lift
                    if movStart.nwash(1,itrial) < movEnd.nwash(1,itrial)
                        % this is good, let's check there is not a lift
                        % between this lift and the grab.
                        liftsBeforeGrab = find(handLift.t0{1,movIdx.nwash(itrial)}{1,1} < movEnd.nwash(1,itrial));
                        if liftsBeforeGrab > 1
                            %this means there are more than one lifts
                            %before the grab, we should choose the one
                            %closest to the grab
                            liftClosest2Grab = handLift.t0{1,movIdx.nwash(itrial)}{1,1}(length(liftsBeforeGrab));
                            movStart.nwash(1,itrial) = liftClosest2Grab;
                        end % if it is one, then that is the only lift
                    else % if the grab happens before the lift, this is an early grab. Find the grab after the lift.
                        % find grab after lift and after digits together
                        firstGrabAfterLift = find(grab.t0{1,movIdx.nwash(itrial)}{1,1} > handLift.t0{1,movIdx.nwash(itrial)}{1,1}(1));
                        % if there is not a grab after lift, either animal
                        % did not finish task or APT did not detect
                        if isempty(firstGrabAfterLift) == 0 %if there is a grab
                        movEnd.nwash(1,itrial) = grab.t0{1,movIdx.nwash(itrial)}{1,1}(firstGrabAfterLift(1));
                        else %if there is no grab after, see if trial is no reach
                           %if reach or grooming, make both zero and move
                           %on
                         %   if trialIdxs.trialScore(base.trialIdxs) == -1
                         % if not reach or grooming, this is an APT error
                         % or add in more conditions
                        end

                    end
                end
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Label grab start as movEnd
                grabStart.nbase = movEnd.nbase;
                grabStart.npert = movEnd.npert;
                grabStart.nwash = movEnd.nwash;
            end % JAABA
            %% Plot Overlay 1D - only reach
            if plot1DOverlayTrajectories == 1
                % Define deminsion to plot (x,y,z)
                dimPlot = 3; % if change add dimPlot input to plotOverlay1D function
                % Plot overlay for baseline
                plotOverlay1D(movStrt.nbase, movEnd.nbase,movIdx.nbase,traj);
                title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' baseline Trials (Dim ' num2str(dimPlot) ')']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayReachBaseline_dim' num2str(dimPlot) '.fig']);
                % Plot overlay for stimulation
                plotOverlay1D(movStrt.npert, movEnd.npert,movIdx.npert,traj);
                title([SUB ' ' EXPER_COND ': Overlay' SCORE ' stim Trials (Dim ' num2str(dimPlot) ')']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayReachStim_dim' num2str(dimPlot) '.fig']);
                % Plot overlay for all baseline, pre and post stimulation
                movStrtBL = [movStrt.nbase movStrt.nwash];
                movEndBL = [movEnd.nbase movEnd.nwash];
                movStrtIdxs = [movIdx.nbase movIdx.nwash]';
                plotOverlay1D(movStrtBL, movEndBL,movStrtIdxs,traj);
                title([SUB ' ' EXPER_COND ': Overlay ' SCORE 'pre and post stim (Dim ' num2str(dimPlot) ')']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayReachAllBL_dim' num2str(dimPlot) '.fig']);
                %% Plot Overlay 1D - whole movement
                % Plot baseline
                plotOverlay1D(movStrt.nbase, movEnd.nbase,movIdx.nbase, traj,'PlotEntireReach',1);
                title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' baseline Trials (Dim ' num2str(dimPlot) ')']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayBaseline_dim' num2str(dimPlot) '.fig']);
                % Plot stimulation
                plotOverlay1D(movStrt.npert, movEnd.npert,movIdx.npert,traj,'PlotEntireReach',1);
                title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' stim Trials (Dim ' num2str(dimPlot) ')']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayStim_dim' num2str(dimPlot) '.fig']);
                % Plot all baseline, pre and post stimulation
                movStrtBL = [movStrt.nbase movStrt.nwash];
                movEndBL = [movEnd.nbase movEnd.nwash];
                movStrtIdxs = [movIdx.nbase movIdx.nwash]';
                plotOverlay1D(movStrtBL, movEndBL,movStrtIdxs,traj,'PlotEntireReach',1);
                title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' pre and post stim (Dim ' num2str(dimPlot) ')']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE 'overlayAllBL_dim' num2str(dimPlot) '.fig']);
            end
            %% Plot Overlay 3D - reach movement
            save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_Overlay3DVariables.mat'],'movStrt', 'movEnd','movIdx', 'traj');
            if plot3DOverlayTrajectories == 1
                % Plot baseline
                plotOverlay3D(movStrt.nbase, movEnd.nbase,movIdx.nbase, traj);
                title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' baseline Trials']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/3DTrajectories/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayReachBaseline.fig']);
                % Plot all baseline, pre and post stimulation
                movStrtBL = [movStrt.nbase movStrt.nwash];
                movEndBL = [movEnd.nbase movEnd.nwash];
                movStrtIdxs = [movIdx.nbase movIdx.nwash]';
                plotOverlay3D(movStrtBL, movEndBL,movStrtIdxs,traj);
                title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' pre and post stim']);
                savefig([BASEPATH 'Figures/' SUB '/Behavior/3DTrajectories/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE 'overlayReachAllBL.fig']);
            end
            %% Get path length for baseline, stim, and washout trials
            % Get path lengths for baseline trials
            [totalPathLength.nbase] = getPathLength(traj, movIdx.nbase,movStrt.nbase, movEnd.nbase);
            % Get path lengths for stimulation trials
            [totalPathLength.npert] = getPathLength(traj, movIdx.npert,movStrt.npert, movEnd.npert);
            % Get path lengths for washout baseline trials
            [totalPathLength.nwash] = getPathLength(traj, movIdx.nwash,movStrt.nwash, movEnd.nwash);
            % Save path lengths in mat file
            save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_PathLengths.mat'],'totalPathLength');
            %% Get digit spread 
            % get pointer digit at first grab for each trial
            % traj(trials x coordinate xyz x frames)
            % find traj points for start of grab
            trajDig1.nbase = traj(movIdx.nbase,:,grabStart.nbase);
            trajDig1.npert = traj(movIdx.npert,:,grabStart.npert);
            trajDig1.nwash = traj(movIdx.nwash,:,grabStart.nwash);
            trajDig2.nbase = trajDigTwo(movIdx.nbase,:,grabStart.nbase);
            trajDig2.npert = trajDigTwo(movIdx.npert,:,grabStart.npert);
            trajDig2.nwash = trajDigTwo(movIdx.nwash,:,grabStart.nwash);

            % baseline trial spread
            for i = 1:length(movIdx.nbase)
                dig1 = trajDig1.nbase(i,:,i);
                dig2 = trajDig2.nbase(i,:,i);
                digitSpread.nbase(i) = sqrt(sum((dig1-dig2) .^ 2));
            end
            % stim trial spread
            for i = 1:length(movIdx.npert)
                dig1 = trajDig1.npert(i,:,i);
                dig2 = trajDig2.npert(i,:,i);
                digitSpread.npert(i) = sqrt(sum((dig1-dig2) .^ 2));
            end
            % wash trial spread
            for i = 1:length(movIdx.nwash)
                dig1 = trajDig1.nwash(i,:,i);
                dig2 = trajDig2.nwash(i,:,i);
                digitSpread.nwash(i) = sqrt(sum((dig1-dig2) .^ 2));
            end    
            save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_DigitSpread.mat'],'digitSpread','movIdx','trialIdxs');
            

        end %session
    end %experimental conditions
end % subject
%% Compare 3D trajectories, harm vs control
Behavioral_Comparison_Conditions; %(animals,BASEPATH, RAWDATA_BASEPATH,exper_conditions,SCORE,USER);
%% Compare accuracy across conditions, harm vs control
Behavioral_Metrics(animals,BASEPATH,exper_conditions);
%% Compare digit spread (pointer and ring), harm vs control
% Calculate distance between digits
% at grab start using traj
Behavioral_DigitSpread; %trajDigSpread, movIdx, trialIdxs);