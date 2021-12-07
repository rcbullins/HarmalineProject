% Neural Comprehensive Script
% PURPOSE
% DEPENDENCIES
%   Run Behavioral Comp Script first (need trajectory mat files)
%   C Drive and google drive have analyzed mat files, code, and figures.
%   D Drive has raw data, videos, execel sheets of recording information.
%   SpikeGLX functions, see function DemoReadSGLXData https://billkarsh.github.io/SpikeGLX/
% OUTPUTS
% HISTORY
%   11.23.2021 Reagan: Adapted from Kevin & Britton
%% Clean workspace
clear;
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
CODE_BRITTON = [RAWDATA_BASEPATH 'Code/britton_code/code/matlab_britton/'];
CODE_SPIKE_GLX = [RAWDATA_BASEPATH 'Code/SpikeGLX_Datafile_Tools/'];
addpath(genpath(CODE_REAGAN));
addpath(genpath(CODE_BRITTON));
addpath(genpath(CODE_SPIKE_GLX));
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
            %% Read lfp data
            lfp_file = [BASEPATH 'Data_Analyzed/' SUB '/Neural/' SUB '_' EXPER_SESSION '_' EXPER_COND '_lfp.mat'];
            if exist(lfp_file) == 0
                % Ask user for binary file
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
                    lfp_chunk_file = [BASEPATH 'Data_Analyzed/' SUB '/Neural/tmp_mat_files/' SUB '_' EXPER_SESSION '_lfpChunk' num2str(ichunk) '.mat'];
                    if exist(lfp_chunk_file) == 0
                        %get waitbar while loading lfp
                        waitbar(fractionVec(ichunk),f,['Loading LFP: ' SUB '-' EXPER_SESSION]);
                        %chunk data into start and stop sample points
                        nSamp_start = nSamp_chunk(1,ichunk);
                        % if on last chunk, go to end
                        if ichunk == size(nSamp_chunk,2)-1
                            nSamp_end = nSamp_chunk(1,ichunk+1);
                        else %if not the end
                            nSamp_end = nSamp_chunk(1,ichunk+1)-1;
                        end
                        % dataArray numChan X samples
                        lfp_chunk = ReadBin(nSamp_start, nSamp_end, meta, ephys_filename,EPHYS_PATH);
                        % For an analog channel: gain correct saved channel ch (1-based for MATLAB).
                        ch = 1;
                        if strcmp(meta.typeThis, 'imec')
                            lfp_chunk = GainCorrectIM(lfp_chunk, [ch], meta);
                        else
                            lfp_chunk = GainCorrectNI(lfp_chunk, [ch], meta);
                        end
                        % downsample lfp from 30000 to 1250
                        ds_factor = length(lfp_chunk)/30000; % lfp once per second
                        ds_factor = ds_factor*1250; % lfp 1250 Hz
                        lfp_downsampled = downsample(lfp_chunk,ds_factor);
                        save([BASEPATH 'Data_Analyzed/' SUB '/Neural/tmp_mat_files/' SUB '_' EXPER_SESSION '_lfpChunk' num2str(ichunk) '.mat'],'lfp_downsampled');
                    else
                        load(lfp_chunk_file);
                    end
                    lfp_all = [lfp_all lfp_downsampled];
                end
                close(f);
                %plot(lfp(ch,:));
                lfp = lfp_all;
                save([BASEPATH 'Data_Analyzed/' SUB '/Neural/' SUB '_' EXPER_SESSION '_' EXPER_COND '_lfp.mat'],'lfp');
                % delete the chunk files once all saved
                for ichunkDel = 1:size(nSamp_chunk,2)-1
                    delete([BASEPATH 'Data_Analyzed/' SUB '/Neural/tmp_mat_files/' SUB '_' EXPER_SESSION '_lfpChunk' num2str(ichunk) '.mat']);
                end
            end
            %% filter lfp - takeout noise (start saving .lf)
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
            if exist([dir_ctx '/../ind.mat']) == 0
                 disp('ind file amiss');
                continue;
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


            %% plot heatmaps of each trial
        
            for n = 1:length(spikes)
                tmp=squeeze(z_norm_lift_ctx(n,:,:));
                subplot(4,6,n);
                hold on
                imagesc(tmp)
                title(strcat('neuron:',num2str(n)));
               yline(length(movStrt.nbase)+1,'--w');
                yline(length(movStrt.nbase)+length(movStrt.npert)+1,'--w');
                xlabel('Time relative to lift')
                xticks([0 500 1000 1500 2000 2500]);
                xticklabels({'-500','0','500','1000','1500','2000','2500'});
            end
            sgtitle([SUB '_' EXPER_SESSION ': ' EXPER_COND ' session firing rates']);

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
             title([SUB '_' EXPER_SESSION ': ' EXPER_COND ' session']);
            s=diag(s);
            VAR_total=sum(s.^2);
            VAF=s.^2/VAR_total;
            %% ISI 
            figure;
            for ineuron = 1:length(spikes)
                thisNeuronSpikes = spikes{ineuron};
                isi=diff(thisNeuronSpikes');
                subplot(4,6,ineuron);
                hist(isi(:),[0:1:2500]);
                xlim([0 1000]);
                title(['Neuron ' num2str(ineuron)]);
            end
            xlabel('Count');
            ylabel('ISI (ms)');
            sgtitle({'Interspike Intervals',[SUB '-' EXPER_SESSION '-' EXPER_COND]});
            %% Autocorr on spikes
            bin_sizes = {100,50,30,25};
            for ibin = 1:length(bin_sizes)
                figure;
                this_bin_size = bin_sizes{ibin};
                for ineuron = 1:length(spikes)
                    thisNeuronSpikes = spikes{ineuron};
                    [N_200,X_200] = hist(reshape(thisNeuronSpikes,1,numel(thisNeuronSpikes)),0:this_bin_size:rec_length*1000);
                    [C,LAGS] = xcov(N_200,'coeff');
                    subplot(4,6,ineuron)
                    plot(LAGS*this_bin_size,C);
                    title(['Neuron ' num2str(ineuron)]);
                end
                sgtitle({['Autocorrelograms: Bin Size ' num2str(this_bin_size) 'ms'],[SUB '-' EXPER_SESSION '-' EXPER_COND]});
            end
            %% Autocorrelogram on lfp
%                 [C,LAGS] = xcov(data(1,:),'coeff');
%                 figure
%                 subplot(2,1,1)
%                 plot(LAGS*.004,C)
%                 xlim([-1 1])
%                 xlabel('Time Lag (s)')
%                 title('Channel 1')
%                 set(gca,'fontsize',32,'fontweight','b')
%% Power spec on lfp
%             sampling_frequency=250; % Set Sampling frequency(line 1)
%             dt=1/sampling_frequency; % Sampling interval (line 2)
%             rec_time=length(data(1,:))./250; % total recording time (line 3)
%             freq_dat_1= fft(data(1,:)); %Calculate the Fourier Transform (line 4)
%             Pow_spec=(2*(dt^2)/rec_time)*abs(freq_dat_1); %Do the calculation for
%             Pow_spec2=Pow_spec(1:(length(data(1,:))/(2))+1); % Only use the
%             df=1/max(rec_time); % Frequency resolution (line 7)
%             fnq= sampling_frequency/2; % Nyquist frequency= half the sampling
%             frequency. (line 8)
%             freq_axis=(0:df:fnq); %Gives you the frequency axis. (line 9)
%             plot(freq_axis,Pow_spec2) % Plot power spectrum (line 10)
%             xlim([0 80]) % Set x-axis limits (line 11)
%             xlabel('Frequency Hz') % (line 12)
%             ylabel ('Power') %(line 13)



        end % experimental session
    end % experimental condition
end % subject




