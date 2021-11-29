% Neural Comprehensive Script
% PURPOSE
% DEPENDENCIES
%   Run Behavioral Comp Script first (need trajectory mat files)
%   C Drive and google drive have analyzed mat files, code, and figures.
%   D Drive has raw data, videos, execel sheets of recording information.
% OUTPUTS
% HISTORY
%   11.23.2021 Reagan: Adapted from Kevin & Britton
%% Clean workspace
clear ALL;
close ALL;
clc;
%% Specify what to plot
score = 'all'; % Options: 1, 0, 2, -1, 'all'
% Code: (1) one grab and success
%       (0) grab and failure
%       (2) multiple reaches and eventual success
%      (-1) no reach attempts
%   ('all') all scores where some attempt was made
%% Add code paths
USER = 'bullinsr'; %'rcbul';
BASEPATH = ['C:/Users/' USER '/OneDrive - University of North Carolina at Chapel Hill/Hantman_Lab/Harmaline_Project/'];
RAWDATA_BASEPATH = 'D:/rbullins/';
CODE_REAGAN = [BASEPATH 'Code/reagan_code/'];
CODE_BRITTON = [RAWDATA_BASEPATH 'Code/britton_code/code'];
addpath(genpath(CODE_REAGAN));
addpath(genpath(CODE_BRITTON));
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
%% For each subject, each condition, each experiment, plot ephys data
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
            %% data sets and trials where magic happens
            % Set file names and directories.
            RAW_EPHYS_FILE = [RAWDATA_BASEPATH 'Data/' SUB 'necab1_Chr2/' EXPER_SESSION 'ephys/' SUB '_' EXPER_SESSION '_1000_g0/' ...
                SUB '_' EXPER_SESSION '_1000_g0_imec0/' SUB '_' EXPER_SESSION '_1000_g0_t0.imec0.ap.bin'];
            RAW_ANALOG_FILE = [RAWDATA_BASEPATH 'Data/' SUB 'necab1_Chr2/' EXPER_SESSION 'ephys/' SUB '_' EXPER_SESSION '_1000_g0/' ...
                SUB '_' EXPER_SESSION '_1000_g0_t0.nidq.meta'];
            %% Set trials
            % Define Trials by baseline, stim, washout
            trials = eval(sprintf('%s_%s_%sTrials',SUB,EXPER_SESSION,EXPER_COND));
            trialIdxs.nbase = trials.nbase;
            trialIdxs.npert = trials.npert;
            trialIdxs.nwash = trials.nwash;
            trialIdxs.trialScore = trials.trialScore;
            %% Parameters
            foo = strfind(RAW_EPHYS_FILE,'/');
            dir_ctx = RAW_EPHYS_FILE(1:(foo(end)-1));

            g_fr = 20; % window size for gaussian smoothing functions
            causal_flag=0;
            window_lift = [-500 2000]; % select data set here
            %% Read metadata (number of physical channels, sampling rate).
            metafile = [RAW_EPHYS_FILE(1:(end-3)) 'meta'];
            meta_text = fileread(metafile);

            %length of recording
            [foo, bar] = regexp(meta_text,'fileTimeSecs=\d');
            rec_length = str2double(meta_text(bar+(0:7))); %length in seconds I believe

            % Read from analog file
            analogMeta_text = fileread(RAW_ANALOG_FILE);
            [foo_sr, bar_sr] = regexp(analogMeta_text,'niSampRate=\d'); %HERE WHAT IS THIS FILE CALLED
            analogSampRate = str2double(analogMeta_text(bar_sr+(0:7)));
            % Set camera sampling rate convert number
            cameraSampRate = 2; %(500 Hz)
            %% Load event indices from the two recording systems and align to Whisper timebase.
            load([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_trajectories.mat'],'traj');
            % Find start frame of movement (uses smoothing and thresholding of velocity)
            % Only for trials with a score of 'score' variable
            [movStrt.nbase] = getMovStartFrame(score,traj,trialIdxs.trialScore,trialIdxs.nbase);
            [movStrt.npert] = getMovStartFrame(score,traj,trialIdxs.trialScore,trialIdxs.npert);
            [movStrt.nwash] = getMovStartFrame(score,traj,trialIdxs.trialScore,trialIdxs.nwash);
            % Find end frame of reach movement (uses smoothing and thresholding of velocity)
            % movIdx is the index of the trial within all trials
            [movEnd.nbase,movIdx.nbase] = getMovEndFrame(score,traj,trialIdxs.trialScore,trialIdxs.nbase, movStrt.nbase);
            [movEnd.npert,movIdx.npert] = getMovEndFrame(score,traj,trialIdxs.trialScore,trialIdxs.npert, movStrt.npert);
            [movEnd.nwash,movIdx.nwash] = getMovEndFrame(score,traj,trialIdxs.trialScore,trialIdxs.nwash, movStrt.nwash);
            % merge all trials together - for now
            startFrames = [movStrt.nbase movStrt.npert movStrt.nwash];
            endFrames = [movEnd.nbase movEnd.npert movEnd.nwash];
            selectIdx = [movIdx.nbase movIdx.npert movIdx.nwash];
            % import ind data
            ind_ctx = load([dir_ctx '/../ind.mat']);
            trial_start_timestamps=get_trial_start(ind_ctx,3)/analogSampRate; %sampling frequency is 16.949 %  IN milliseconds
            % Select indexes for interest
            trial_start_timestamps_idxSelect = trial_start_timestamps(selectIdx);
            % Get movement onset time (start)
            trial_start_mov_timestamps = trial_start_timestamps_idxSelect + (startFrames.*cameraSampRate);

            %% Get spikes.
            % Cortical spikes.
            [spk, mua, in] = get_st_ks2(dir_ctx);
            for j = 1:length(spk)
                st_ctx{j} = double(spk{j})/30; %divide by 30 because neural signals sampled at 30kHz
            end
            %% Peri-lift firing rates
            clear mu_ctx std_ctx

            % Get cortical firing rates.
            [rates,t_rates]=convolve_time_stamps(st_ctx,g_fr,rec_length*1000,causal_flag);

            %% cut up time-series into each trial
            r_lift_ctx=[]; %neuron x trial x timepoints


            for t=1:length(trial_start_mov_timestamps)
                wind_min=round(trial_start_mov_timestamps(t)+window_lift(1));
                % Want to go back in space, but if trial starts before that
                % window, just take out trial (R fix)
                if wind_min < 0 
                    wind_min = 1;
                    r_lift_ctx(:,t,:) = NaN(size(rates,1),(abs(window_lift(1))+window_lift(2))+1);
                    continue;
                end
                wind_max=round(trial_start_mov_timestamps(t)+window_lift(2));
                r_lift_ctx(:,t,:)=rates(:,wind_min:wind_max);
            end
            % Take out trials if not enough time points before it
                %%%%HERE TAKE OUT NAN ROWS%%%%%%%

            %% z-score firing rates
            % Get z-scores, WITH SOFT NORMALIZATION?
            i4z = 1:400;
            c_softnorm = .5;
            c_soft_z = 5;
            z_norm_lift_ctx=[];
            for i = 1:size(r_lift_ctx,1)
                dd = squeeze(r_lift_ctx(i,:,i4z));
                dd = reshape(dd,1,numel(dd));
                mu_ctx(i) = nanmean(dd);
                std_ctx(i) = nanstd(dd);

                %lift aligned
                z_norm_lift_ctx(i,:,:) = (squeeze(r_lift_ctx(i,:,:))-mu_ctx(i))./(c_soft_z + std_ctx(i));
            end


            %% select trials of interest
            % Baseline
            z_norm_lift_ctx=z_norm_lift_ctx(:,trialIdxs,:);


            %% plot heatmaps of each trial
            for n=1:2:10
                tmp=squeeze(z_norm_lift_ctx(n,:,:));
                figure
                hold on
                imagesc(tmp)
                title(strcat('neuron:',num2str(n)))
                xlabel('Time relative to lift')
            end

            %% perform PCA analysis
            n_pc=3;
            %trial average
            z_norm_tavg=squeeze(mean(z_norm_lift_ctx,2));

            % subtract off mean activity
            %average across time
            mean_z_norm_tavg=squeeze(mean(z_norm_tavg,2)); %neurons x 1
            %expand it to be neurons x time
            ntimepts=size(z_norm_tavg,2);
            mean_z_norm_tavg=repmat(mean_z_norm_tavg,[1,ntimepts]);

            z_norm_tavg=z_norm_tavg-mean_z_norm_tavg;

            % Get PCA coefficients from trial averaged
            [u,s,v]=svd(z_norm_tavg); % all trials and timepoints
            PC=u(:,1:3).';

            proj_data=PC*z_norm_tavg;


            figure
            hold on
            plot3(proj_data(1,:),proj_data(2,:),proj_data(3,:))
            xlabel('PC 1')
            ylabel('PC 2')
            zlabel('PC 3')

            s=diag(s);
            VAR_total=sum(s.^2);
            VAF=s.^2/VAR_total;
        end % experimental session
    end % experimental condition
end % subject




