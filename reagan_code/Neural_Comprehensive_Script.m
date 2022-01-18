% Neural Comprehensive Script
% PURPOSE
% DEPENDENCIES
%   Run APT > Run Behavioral_Tracking > Run JAABA GUI > Output predictions
%   Run Behavioral Comp Script first (need trajectory mat files)
%   C Drive and google drive have analyzed mat files, code, and figures.
%   D Drive has raw data, videos, execel sheets of recording information.
%   SpikeGLX functions, see function DemoReadSGLXData https://billkarsh.github.io/SpikeGLX/
% OUTPUTS
% HISTORY
%   11.23.2021 Reagan: Some sections adapted from Kevin & Britton
%% Clean workspace
clear;
close ALL;
clc;
%% Run specifics
useJAABA     = 0; %use JAABA output to find lift times (1) or APT thresholding (0)
runSpikeInfo = 0;
runLFPInfo   = 1;
runHeatMaps  = 0;
runPCA       = 0;
runDecoder   = 0;
runISI       = 0;
runAutoCorr  = 0;
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
JAABA_OUTPUT = [RAWDATA_BASEPATH 'JAABA_behaviors/'];
CODE_REAGAN = [BASEPATH 'Code/'];
CODE_BRITTON = [RAWDATA_BASEPATH 'Code/britton_code/code/matlab_britton/'];
CODE_SPIKE_GLX = [RAWDATA_BASEPATH 'Code/SpikeGLX_Datafile_Tools/'];
CODE_LFP_ANALYSIS = [RAWDATA_BASEPATH 'Code/PlaceInhibitionCode/'];
addpath('D:\rbullins\Code\decoder');
PROBES = 'D:\rbullins\Code\process_events\probes';
addpath(genpath(PROBES));
addpath(genpath(CODE_LFP_ANALYSIS));
addpath(genpath(CODE_REAGAN));
addpath(genpath(CODE_BRITTON));
addpath(genpath(CODE_SPIKE_GLX));
% Add paths to data inventory script (edit to specify sessions to run over)
Directory_Animals;
%% Experimental conditions
exper_conditions = {'control';'harm'};

colorMap(1,:) =[70/255 130/255 180/255];
colorMap(2,:) = [178/255 34/255 34/255];
%% Classifiers
classifiers = {'handLift', 'digitsTogether', 'grab', 'supinate','atMouth'};
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
%% Set Graph Defaults - I like arial font :)
SetGraphDefaults;
%% For each subject, each condition, each experiment, plot ephys data
% Loop through subjects
for isub = 1%:length(animals)
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
            RAW_IND_BIN = [RAWDATA_BASEPATH 'Data/' SUB 'necab1_Chr2/' EXPER_SESSION 'ephys/' SUB '_' EXPER_SESSION '_1000_g0/'...
                SUB '_' EXPER_SESSION '_1000_g0_t0.nidq.bin'];
            if strcmp(EXPER_SESSION,'20210819')
                RAW_ANALOG_FILE = 'D:/rbullins/Data/M341necab1_Chr2/20210819ephys/M341_20210819_1100_g0/M341_20210819_1100_g0_t0.nidq.meta';
                RAW_EPHYS_FILE = 'D:/rbullins/Data/M341necab1_Chr2/20210819ephys/M341_20210819_1100_g0/M341_20210819_1100_g0_imec0/M341_20210819_1100_g0_t0.imec0.ap.bin';
                RAW_IND_BIN = 'D:/rbullins/Data/M341necab1_Chr2/20210819ephys/M341_20210819_1100_g0/M341_20210819_1100_g0_t0.nidq.bin';
            end
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
            analogSampRate = str2double(analogMeta_text(bar_sr+(0:7))); %i think
            % Set camera sampling rate convert number
            cameraSampRate = 2; %(500 Hz)
            if useJAABA
                %% Load classifier data (hand lift, grab, etc)
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
            end

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
            for iclassifier = 1:length(classifiers)
                CLASS = classifiers{iclassifier};
                if useJAABA
                    startFrames = zeros(1,length(selectIdx));
                    for itrial = 1:length(selectIdx) %1:numTrials
                        thisClass = eval(sprintf('%s',CLASS));
                        if isempty(thisClass.t0{1,selectIdx(itrial)}{1,1}) == 0
                            startFrames(1,itrial) = thisClass.t0{1,selectIdx(itrial)}{1,1}(1);
                        else
                            startFrames(1,itrial) = 0;
                        end
                    end
                end
                % import ind data
                if exist([dir_ctx '/../ind.mat'])==0
                    %make ind file
                    f_in = RAW_IND_BIN;
                    foo = strfind(f_in,'/');
                    dir_out = f_in(1:(foo(end)-1));
                    thresh = [.001 .002 .001 .002 .002];

                    % ( fname,nchan,chan2scan,thresh,invert )
                    % ind = get_event_ind(f_in,67,65:67,thresh,0);
                    cd(PROBES)
                    ind = get_event_ind(f_in,5,[1:5],thresh,0);
                    cd('../')
                    % i_first_cam=ind{5};

                    save([dir_out '/ind.mat'],'ind');
                end
                ind_ctx = load([dir_ctx '/../ind.mat']);
                trial_start_timestamps=get_trial_start(ind_ctx,3)/analogSampRate; %sampling frequency is 16.949 %  IN seconds!!!
                % Select indexes for interest
                trial_start_timestamps_idxSelect = trial_start_timestamps(selectIdx);
                % Get movement onset time (start)
                trial_start_mov_timestamps = trial_start_timestamps_idxSelect + ((startFrames.*cameraSampRate)/1000); % In seconds
                trial_start_mov_timestamps = trial_start_mov_timestamps.*1000; % convert to miliseconds

                %% Get spikes.
                % Load spikes from kilosort output
                [spk, mua, in] = get_st_ks2(dir_ctx);
                for j = 1:length(spk)
                    spikes{j} = double(spk{j})/30; %divide by 30 because neural signals sampled at 30kHz, no in ms
                end

                %% Peri-lift firing rates
                clear mu_ctx std_ctx

                % Get cortical firing rates.
                [rates,t_rates]=convolve_time_stamps(spikes,g_fr,rec_length*1000,causal_flag); %in ms

                %% cut up time-series into each trial
                r_lift_ctx=[]; %neuron x trial x timepoints
                takeOutIdx = [];
                % for each trial, take time before trial start of window 1 and
                % take samples until window 2
                for t=1:length(trial_start_mov_timestamps)
                    wind_min=round(trial_start_mov_timestamps(t)+window_lift(1));
                    % Want to go back in space, but if trial starts before that
                    % window, just take out trial (R fix)
                    if wind_min < 0
                        wind_min = 1;
                        r_lift_ctx(:,t,:) = NaN(size(rates,1),(abs(window_lift(1))+window_lift(2))+1);
                        takeOutIdx = [takeOutIdx t];
                        disp(['taking out trial index ' num2str(t)]);
                        continue;
                    end
                    wind_max=round(trial_start_mov_timestamps(t)+window_lift(2));
                    r_lift_ctx(:,t,:)=rates(:,wind_min:wind_max);
                end
                % Take out trials if not enough time points before it
                r_lift_ctx(:,takeOutIdx,:) = [];

                %% z-score firing rates
                % Get z-scores, WITH SOFT NORMALIZATION?
                i4z = 1:400;
                c_softnorm = .5;
                c_soft_z = 5;
                z_norm_lift_ctx=[];
                for i = 1:length(spikes)
                    dd = squeeze(r_lift_ctx(i,:,i4z));
                    dd = reshape(dd,1,numel(dd));
                    mu_ctx(i) = nanmean(dd);
                    std_ctx(i) = nanstd(dd);

                    %lift aligned
                    z_norm_lift_ctx(i,:,:) = (squeeze(r_lift_ctx(i,:,:))-mu_ctx(i))./(c_soft_z + std_ctx(i));
                end


                %% select trials of interest
                % Define Trials by baseline, stim, washout
                trialIdxs = eval(sprintf('%s_%s_%sTrials',SUB,EXPER_SESSION,EXPER_COND));
                %selectIdx
                % Baseline
                z_norm_lift_ctx=z_norm_lift_ctx(:,:,:); % second dim is index of wanted trials
                % for subplots in plotting
                totalNumCells = length(spikes);
                columnNumCells = 6;
                rowsNumCells = ceil(totalNumCells/columnNumCells);


                if runHeatMaps
                    %% plot heatmaps of each trial
                    % for size of subplot
                    for n = 1:length(spikes)
                        tmp=squeeze(z_norm_lift_ctx(n,:,:));
                        subplot(rowsNumCells,columnNumCells,n);
                        hold on
                        imagesc(tmp)
                        sgtitle(strcat('neuron:',num2str(n)));
                        yline(length(movStrt.nbase)+1,'--w');
                        yline(length(movStrt.nbase)+length(movStrt.npert)+1,'--w');
                        xlabel('Time')
                        xticks([0 500 1000 1500 2000 2500]);
                        xticklabels({'-500','0','500','1000','1500','2000','2500'});
                    end
                    sgtitle({[SUB ' ' EXPER_SESSION ': ' EXPER_COND ' session firing rates'], ['Aligned to ' CLASS]});
                    savefig([BASEPATH 'Figures/' SUB '/Neural/' SUB '_' EXPER_SESSION '_' EXPER_COND '_FR_' CLASS '.fig']);
                    %% Plot heatmap of all neurons
                    figure;
                    heatMapAll = zeros(length(spikes),size(z_norm_lift_ctx,3));
                    for ineuron = 1:length(spikes)
                        tmp = squeeze(z_norm_lift_ctx(ineuron,:,:));
                        thisNeuron = mean(tmp);
                        heatMapAll(ineuron,:) = thisNeuron;
                    end
                    imagesc(heatMapAll);
                    xlabel('Time')
                    ylabel('Neuron')
                    xticks([0 500 1000 1500 2000 2500]);
                    xticklabels({'-500','0','500','1000','1500','2000','2500'});
                    sgtitle({[SUB ' ' EXPER_SESSION ': ' EXPER_COND ' session firing rates'], ['Aligned to ' CLASS]});
                    savefig([BASEPATH 'Figures/' SUB '/Neural/' SUB '_' EXPER_SESSION '_' EXPER_COND '_FR_Avg_' CLASS '.fig']);
                end
                if runPCA
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


                    figure;
                    hold on;
                    plot3(proj_data(1,:),proj_data(2,:),proj_data(3,:));
                    xlabel('PC 1');
                    ylabel('PC 2');
                    zlabel('PC 3');
                    sgtitle({[SUB ' ' EXPER_SESSION ': ' EXPER_COND ' session'], ['Aligned to ' CLASS]});
                    s=diag(s);
                    VAR_total=sum(s.^2);
                    VAF=s.^2/VAR_total;
                    savefig([BASEPATH 'Figures/' SUB '/Neural/' SUB '_' EXPER_SESSION '_' EXPER_COND '_PCA_' CLASS '.fig']);
                end %PCA
            end % classifiers
            if runISI
                %% ISI
                figure;
                peakISI.ms = zeros(1,length(spikes));

                for ineuron = 1:length(spikes)
                    thisNeuronSpikes = spikes{ineuron};
                    isi=diff(thisNeuronSpikes');
                    subplot(rowsNumCells,columnNumCells,ineuron);
                    hgram = histogram(isi(:),[0:1:max(isi)]);%2500]);
                    xlim([0 1000]);% up to one second interval
                    title(['Neuron ' num2str(ineuron)]);
                    maxPeak = find(hgram.Values==max(hgram.Values))+1;
                    if length(maxPeak) > 1
                        maxPeak = round(mean(maxPeak));
                    end
                    peakISI.ms(1,ineuron)= hgram.BinEdges(maxPeak); %in ms

                end
                xlabel('Count');
                ylabel('ISI (ms)');
                sgtitle({'Interspike Intervals',[SUB '-' EXPER_SESSION '-' EXPER_COND]});

                savefig([BASEPATH 'Figures/' SUB '/Neural/' SUB '_' EXPER_SESSION '_' EXPER_COND '_ISI.fig']);

                % Second ISI
                figure;
                hist(peakISI.ms,(0:1:max(peakISI.ms)));
                xlabel('ISI (ms)');
                ylabel('Count');
                box off;
                sgtitle({'Interspike Intervals Peaks',[SUB '-' EXPER_SESSION '-' EXPER_COND]});

                savefig([BASEPATH 'Figures/' SUB '/Neural/' SUB '_' EXPER_SESSION '_' EXPER_COND '_ISI_Peaks.fig']);
            end
            if runAutoCorr
                %% Autocorr on spikes
                bin_sizes = {100,50,30,25};
                for ibin = 1:length(bin_sizes)
                    figure;
                    this_bin_size = bin_sizes{ibin};
                    for ineuron = 1:length(spikes)
                        thisNeuronSpikes = spikes{ineuron};
                        [N_200,X_200] = hist(reshape(thisNeuronSpikes,1,numel(thisNeuronSpikes)),0:this_bin_size:rec_length*1000);
                        [C,LAGS] = xcov(N_200,'coeff');
                        subplot(rowsNumCells,columnNumCells,ineuron)
                        plot(LAGS*this_bin_size,C);
                        title(['Neuron ' num2str(ineuron)]);
                    end
                    sgtitle({['Autocorrelograms: Bin Size ' num2str(this_bin_size) 'ms'],[SUB '-' EXPER_SESSION '-' EXPER_COND]});

                    savefig([BASEPATH 'Figures/' SUB '/Neural/' SUB '_' EXPER_SESSION '_' EXPER_COND '_AutoCorrSpikes.fig']);
                end
            end
            if runLFPInfo
                % Sampling Frequency
                Fs = 1250;
                %% Read lfp data
                lfp_file = [BASEPATH 'Data_Analyzed/' SUB '/Neural/' SUB '_' EXPER_SESSION '_' EXPER_COND '_lfp.mat'];
                lfp_filtered_file = [BASEPATH 'Data_Analyzed/' SUB '/Neural/' SUB '_' EXPER_SESSION '_' EXPER_COND '_lfpFilt.mat'];
                if exist(lfp_filtered_file)==0
                    if exist(lfp_file)==0
                        lfp = getLFPfromBin(RAW_EPHYS_FILE,BASEPATH, rec_length, SUB, EXPER_SESSION);
                        save([BASEPATH 'Data_Analyzed/' SUB '/Neural/' SUB '_' EXPER_SESSION '_' EXPER_COND '_lfp.mat'],'lfp', '-v7.3');
                    else
                        load(lfp_file,'lfp');
                    end
                    % filter lfp - takeout noise
                    d = designfilt('bandstopiir','FilterOrder',2, ...
                        'HalfPowerFrequency1',58,'HalfPowerFrequency2',62, ...
                        'DesignMethod','butter','SampleRate',Fs);
                    lfp = lfp(:,1:4500000);
                    butterLoop = filtfilt(d,lfp);
                    plot(lfp(1:Fs),'r')
                    hold on
                    plot(butterLoop(1:Fs),'b')
                    lfp_filtered = lowpass(lfp,80,Fs);
                    save([BASEPATH 'Data_Analyzed/' SUB '/Neural/' SUB '_' EXPER_SESSION '_' EXPER_COND '_lfpFilt.mat'],'lfp_filtered', '-v7.3');
                else
                    load(lfp_filtered_file,'lfp_filtered');
                end
                %% Autocorrelogram on lfp
                chan2Plot = 1;
                data = lfp_filtered;
                %data = bandpass(lfp_filtered,[58 62],1250);
                [C,LAGS] = xcov(data(chan2Plot,:),'coeff');
                figure
                subplot(2,1,1)
                plot(LAGS*.004,C)
                xlim([-1 1])
                xlabel('Time Lag (s)')
                sgtitle({['Autocorrelogram Channel: ' num2str(chan2Plot)], [SUB '-' EXPER_SESSION '-' EXPER_COND]})
                set(gca,'fontsize',32,'fontweight','b')

                savefig([BASEPATH 'Figures/' SUB '/Neural/' SUB '_' EXPER_SESSION '_' EXPER_COND '_AutoCorrLFP.fig']);
                %% Power spec on lfp
                dt=1/Fs; % Sampling interval
                rec_time= length(data(chan2Plot,:))./Fs; % total recording time
                freq_dat_1= fft(data(chan2Plot,:)); %Calculate the Fourier Transform
                Pow_spec=(2*(dt^2)/rec_time)*abs(freq_dat_1); %Do the calculation for
                Pow_spec2=Pow_spec(1:(length(data(chan2Plot,:))/(2))+1); % Only use the
                df=1/max(rec_time); % Frequency resolution
                fnq= Fs/2; % Nyquist frequency= half the sampling frequency.
                freq_axis=(0:df:fnq); %Gives you the frequency axis.
                plot(freq_axis,Pow_spec2,'Color',colorMap(iexper,:)); % Plot power spectrum
                xlim([0 50]); % Set x-axis limits
                xlabel('Frequency Hz');
                ylabel ('Power');
                sgtitle({['LFP PowerSpec Channel: ' num2str(chan2Plot)], [SUB '-' EXPER_SESSION '-' EXPER_COND]});

                savefig([BASEPATH 'Figures/' SUB '/Neural/' SUB '_' EXPER_SESSION '_' EXPER_COND '_PowerspecLFP.fig']);
                if strcmp(EXPER_COND, 'control')
                    freq_axis_control = freq_axis;
                    Pow_spec2_control = Pow_spec2;
                    lfp_control_example = data(chan2Plot,75000:75000+1250);
                elseif strcmp(EXPER_COND,'harm')
                    freq_axis_harm = freq_axis;
                    Pow_spec2_harm = Pow_spec2;
                    lfp_harm_example = data(chan2Plot,75000:75000+1250);
                end
                %% LFP legit way
                % load lfp
                % make .data and .timestamps lfp
                lfpStruct.data = lfp_filtered;
                lfpStruct.data = bandpass(lfpStruct.data,[58 62],1250);
                totalSecRec = length(lfpStruct.data)/1250;
                lfpStruct.timestamps = (0:1/Fs:totalSecRec-(1/Fs));
%                 % Get power of data
%                 doSplitLFP = 0;
%                 if ~doSplitLFP
%                     [powS1] = getPowerSpectrum(lfpStruct, 'doIRASA', false,'doPlot', false);
%                 elseif doSplitLFP
%                     minutes = 5;
%                     timeMin = 1250*60*minutes;
%                     % Sleep 1
%                     [powS1] = makePowerspec_AvgChunk_Vector(basePath, timeMin, lfp);
%                 end
                % Plot Power
                                figure;
                                plot(freq_axis,Pow_spec2);
                                legend({'Control'});
                                legend boxoff;
                                title('Powerspectrum');
                                xlabel('Frequency (Hz)');
                                ylabel('Power (mV)');
                                % Chunk lfp and Take out the fractals
                                % take out fractals
%                                 if ~doSplitLFP
                                    time30min = 30*(60*1250);
                                    lfp_channel = 1;
                                    specS1 = amri_sig_fractal(lfpStruct.data(lfp_channel,1:time30min), 1250,'detrend',1,'frange', [1 150]);
%                                 elseif doSplitLFP
%                                     %make movmean 1
%                                     movmean_win = 1;
%                                     [specS1_freq, specS1_osci] = chunkLFP_takeOutFractals(lfp_VR);
%                                     specS1.freq = mean(specS1_freq);
%                                     specS1.osci = mean(specS1_osci);
%                                 end
                                figure;
                                movmean_win = 1000;
                                plot(specS1.freq,movmean(specS1.osci,movmean_win));
                                hold on ;
                                title('Powerspectrum for Sleep Sessions:No fractals');
                                xlabel('Frequency (Hz)');
                                ylabel('Power (mV)');
                                % Make IRASA struct IRASA
                                IRASA.specS1 = specS1;
            end
         
        end % experimental session
    end % experimental condition
    % Plot LFP harm vs control for this animal
    if runLFPInfo
        figure;
        plot(freq_axis_control,Pow_spec2_control,'Color',colorMap(1,:));
        hold on;
        plot(freq_axis_harm,Pow_spec2_harm,'Color',colorMap(2,:));
        legend('Control','Harmaline');
        legend boxoff;
        xlim([0 50]); % Set x-axis limits
        ylim([0 1*10^-7]);
        xlabel('Frequency Hz');
        ylabel ('Power');
        sgtitle({['LFP PowerSpec Channel: ' num2str(chan2Plot)], [SUB '-' EXPER_SESSION]});

        savefig([BASEPATH 'Figures/' SUB '/Neural/' SUB '_' EXPER_SESSION '_' EXPER_COND '_Powerspec_Comp.fig']);
        figure;
        plot(lfp_control_example,'Color',colorMap(1,:));
        hold on;
        plot(lfp_harm_example,'Color',colorMap(2,:));
        legend('Control','Harmaline');
        legend boxoff;
        xlabel('Time (ms)')
        ylabel('mV')
        xticks([0 312.5 625 937.5 1250])
        xticklabels({'0','250','500','750','1000'})
        sgtitle({['LFP Trace Channel: ' num2str(chan2Plot)], [SUB '-' EXPER_SESSION]});

        savefig([BASEPATH 'Figures/' SUB '/Neural/' SUB '_' EXPER_SESSION '_' EXPER_COND '_LFPTraceExample.fig']);

    end
end % subject


