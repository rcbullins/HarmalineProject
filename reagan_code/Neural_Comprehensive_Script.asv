% Neural Comprehensive Script
% PURPOSE
% DEPENDENCIES
% OUTPUTS
% HISTORY
%   11.23.2021 Reagan Adapted from Kevin and Britton
%% Clean workspace
clear ALL;
close ALL;
clc;
%% Specify what to plot
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
USER = 'bullinsr'; %'rcbul';
BASEPATH = ['C:/Users/' USER '/OneDrive - University of North Carolina at Chapel Hill/Hantman_Lab/Harmaline_Project/'];
CODE_REAGAN = [BASEPATH 'Code/reagan_code/'];
CODE_BRITTON = [BASEPATH 'Code/britton_code/code'];
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
        BEHAV_FIG = [BASEPATH 'Figures\' SUB '\Behavior\'];
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
RAW_EPHYS_FILE = [BASEPATH 'Data/' SUB 'necab1_Chr2/' EXPER_SESSION 'ephys/' SUB '_' EXPER_SESSION '_1000_g0/' ...
                 SUB '_' EXPER_SESSION '_1000_g0_imec0/' SUB '_' EXPER_SESSION '_1000_g0_t0.imec0.ap.bin'];
%% Set trials
% Define Trials by baseline, stim, washout
trials = eval(sprintf('%s_%s_%sTrials',SUB,EXPER_SESSION,EXPER_COND));
trialIdxs.base = trials.nbase;
trialIdxs.pert = trials.npert;
trialIdxs.wash = trials.nwash;
%% Parameters
foo = strfind(RAW_EPHYS_FILE,'\');
dir_ctx = RAW_EPHYS_FILE(1:(foo(end)-1));

g_fr = 20; % window size for gaussian smoothing functions
causal_flag=0;
window_lift = [-1000 2000]; % select data set here
%% Read metadata (number of physical channels, sampling rate).
metafile = [RAW_EPHYS_FILE(1:(end-3)) 'meta'];
meta_text = fileread(metafile);

%length of recording
[foo bar] = regexp(meta_text,'fileTimeSecs=\d');
fil_len = str2double(meta_text(bar+(0:7)));
%%%%%% WORKING HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HERE - Rayg figure out how to load sampling freq, and soft code line
% below


%% Load event indices from the two recording systems and align to Whisper timebase.
% Load start and end frame information from behavioral script

% import ind data
ind_ctx = load([dir_ctx '/../ind.mat']);
trial_start_timestamps=get_trial_start(ind_ctx,3)/16.949; %sampling frequency is 16.949
% trial_start = trial_start_timestamps + startFrame


%% Get spikes.
% Cortical spikes.

[spk mua in] = get_st_ks2(dir_ctx);
for j = 1:length(spk)
    st_ctx{j} = double(spk{j})/30; %divide by 30 because neural signals sampled at 30kHz
end



%% Peri-lift firing rates 

clear mu_ctx std_ctx

% Get cortical firing rates.
[rates,t_rates]=convolve_time_stamps(st_ctx,g_fr,fil_len*1000,causal_flag);

%% cut up time-series into each trial
r_lift_ctx=[]; %neuron x trial x timepoints


for t=1:length(trial_start_timestamps)
    wind_min=round(trial_start_timestamps(t)+window_lift(1));
    wind_max=round(trial_start_timestamps(t)+window_lift(2));
    r_lift_ctx(:,t,:)=rates(:,wind_min:wind_max);
end



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




