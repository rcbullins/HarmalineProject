% Testing Decoder
%% Variables
window_lift = [1 5000]; % select data set here
%% Set Paths and Experimental Session
USER = 'bullinsr';
BASEPATH = ['C:/Users/' USER '/OneDrive - University of North Carolina at Chapel Hill/Hantman_Lab/Harmaline_Project/'];
RAWDATA_BASEPATH = 'D:/rbullins/';
JAABA_OUTPUT = [RAWDATA_BASEPATH 'JAABA_behaviors/'];
Directory_Animals;
SetGraphDefaults;
score = [1 0 2];
classifiers = {'handLift'};% 'digitsTogether', 'grab', 'supinate','atMouth'};
CLASS = classifiers{1};
useJAABA = 0;
%% Experimental conditions
exper_conditions = {'control';'harm'};
%% For each subject, each condition, each experiment, plot ephys data
% Loop through subjects
for isub = 2%:length(animals)
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

            trialIdxs = eval(sprintf('%s_%s_%sTrials',SUB,EXPER_SESSION,EXPER_COND));

            RAW_EPHYS_FILE = [RAWDATA_BASEPATH 'Data/' SUB 'necab1_Chr2/' EXPER_SESSION 'ephys/' SUB '_' EXPER_SESSION '_1000_g0/' ...
                SUB '_' EXPER_SESSION '_1000_g0_imec0/' SUB '_' EXPER_SESSION '_1000_g0_t0.imec0.ap.bin'];
            RAW_ANALOG_FILE = [RAWDATA_BASEPATH 'Data/' SUB 'necab1_Chr2/' EXPER_SESSION 'ephys/' SUB '_' EXPER_SESSION '_1000_g0/' ...
                SUB '_' EXPER_SESSION '_1000_g0_t0.nidq.meta'];
            RAW_IND_BIN = [RAWDATA_BASEPATH 'Data/' SUB 'necab1_Chr2/' EXPER_SESSION 'ephys/' SUB '_' EXPER_SESSION '_1000_g0/'...
                SUB '_' EXPER_SESSION '_1000_g0_t0.nidq.bin'];

               if strcmp(SUB,'M341')==1
                   if strcmp(EXPER_SESSION,'20210819')==1
 RAW_EPHYS_FILE = [RAWDATA_BASEPATH 'Data/' SUB 'necab1_Chr2/' EXPER_SESSION 'ephys/' SUB '_' EXPER_SESSION '_1100_g0/' ...
                SUB '_' EXPER_SESSION '_1100_g0_imec0/' SUB '_' EXPER_SESSION '_1100_g0_t0.imec0.ap.bin'];
            RAW_ANALOG_FILE = [RAWDATA_BASEPATH 'Data/' SUB 'necab1_Chr2/' EXPER_SESSION 'ephys/' SUB '_' EXPER_SESSION '_1100_g0/' ...
                SUB '_' EXPER_SESSION '_1100_g0_t0.nidq.meta'];
            RAW_IND_BIN = [RAWDATA_BASEPATH 'Data/' SUB 'necab1_Chr2/' EXPER_SESSION 'ephys/' SUB '_' EXPER_SESSION '_1100_g0/'...
                SUB '_' EXPER_SESSION '_1100_g0_t0.nidq.bin'];
                   end
               end
            %% Load Trajectory
            trajFile = ([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_trajectories.mat']);
            if exist(trajFile,'file')==2
                load(trajFile);
            else
                nframes=2500; %variable! %number of frames to analyze. Going to depend on how many frames in movie
                digit2trk = 1;
                [traj, ~] = get_traj_3D_v2(TRK,CODE_CALIB,nframes,digit2trk);
                save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_trajectories.mat'],'traj');
            end

            %% Load event indices
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
            %             selectIdx = movIdx.nbase;%[movIdx.nbase movIdx.npert movIdx.nwash];
            selectIdx_ideal_bl = trialIdxs.nbase(find(trialIdxs.trialScore(trialIdxs.nbase) == 1));
            selectIdx_eventualSuccess_bl = trialIdxs.nbase(find(trialIdxs.trialScore(trialIdxs.nbase) == 2));
            selectIdx_noSuccess_bl = trialIdxs.nbase(find(trialIdxs.trialScore(trialIdxs.nbase) == 0));
            selectIdx_bl = [selectIdx_ideal_bl]; % selectIdx_eventualSuccess selectIdx_noSuccess];
            selectIdx_ideal_wash = trialIdxs.nwash(find(trialIdxs.trialScore(trialIdxs.nwash) == 1));
            selectIdx_eventualSuccess_wash = trialIdxs.nwash(find(trialIdxs.trialScore(trialIdxs.nwash) == 2));
            selectIdx_noSuccess_wash = trialIdxs.nwash(find(trialIdxs.trialScore(trialIdxs.nwash) == 0));
            selectIdx_wash = [selectIdx_ideal_wash];
            selectIdx = [selectIdx_bl selectIdx_wash];

           % Start frames via APT
           start_ideal_bl = movStrt.nbase(find(trialIdxs.trialScore(movIdx.nbase) == 1));
            start_ideal_wash = movStrt.nwash(find(trialIdxs.trialScore(movIdx.nwash) == 1));
            startFrames_APT = [start_ideal_bl start_ideal_wash];
            %% Get hand lift times from JAABA
if useJAABA
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
            handLift.t0       = {};
            % Initiate cells to store stop frames in
            handLift.t1       = {};
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
                handLift_scores       = load([TRIAL_FOLDER 'scores_LiftHand.mat'],'allScores');
                % save in mega cell array start times
                handLift.t0{itrial}       = handLift_scores.allScores.t0s;
                % save in mega cell array stop times
                handLift.t1{itrial}       = handLift_scores.allScores.t1s;
            end
            %% Get start frames of movement hand lift
            startFrames_handLift = zeros(1,length(selectIdx));

            for itrial = 1:length(selectIdx) %1:numTrials
                thisClass = eval(sprintf('%s',CLASS));
                if isempty(thisClass.t0{1,selectIdx(itrial)}{1,1}) == 0
                    startFrames_handLift(1,itrial) = thisClass.t0{1,selectIdx(itrial)}{1,1}(1);
                else
                    startFrames_handLift(1,itrial) = 0;
                end
            end
end
            %% Params
            foo = strfind(RAW_EPHYS_FILE,'/');
            dir_ctx = RAW_EPHYS_FILE(1:(foo(end)-1));

            g_fr = 20; % window size for gaussian smoothing functions
            causal_flag=0;
            % Read from analog file
         
            analogMeta_text = fileread(RAW_ANALOG_FILE);
            [foo_sr, bar_sr] = regexp(analogMeta_text,'niSampRate=\d'); %HERE WHAT IS THIS FILE CALLED
            analogSampRate = str2double(analogMeta_text(bar_sr+(0:7)))/1000;
            % Set camera sampling rate convert number
            cameraSampRate = 2; %(500 Hz)
            metafile = [RAW_EPHYS_FILE(1:(end-3)) 'meta'];
            meta_text = fileread(metafile);

            %length of recording
            [foo, bar] = regexp(meta_text,'fileTimeSecs=\d');
            rec_length = str2double(meta_text(bar+(0:7))); %length in seconds I believe

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
            %% Get trial start time and timestamp of movement initiation
            trial_start_timestamps=get_trial_start(ind_ctx,3)/analogSampRate; %sampling frequency is 16.949 %  IN seconds!!!
            % Select indexes for interest - trial start times in recording
            trial_start_timestamps_idxSelect = trial_start_timestamps(selectIdx);

            %% Set start times to get neural data from
            startTimes = trial_start_timestamps_idxSelect;

            %% Get spikes.
            % Load spikes from kilosort output
            [spk, mua, in] = get_st_ks2(dir_ctx);
            for j = 1:length(spk)
                spikes{j} = double(spk{j})/30; %divide by 30 because neural signals sampled at 30kHz, now in ms
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
            for t=1:length(startTimes)
                wind_min=round(startTimes(t)+window_lift(1));
                % Want to go back in space, but if trial starts before that
                % window, just take out trial (R fix)
                if wind_min < 0
                    wind_min = 1;
                    r_lift_ctx(:,t,:) = NaN(size(rates,1),(abs(window_lift(1))+window_lift(2))+1);
                    takeOutIdx = [takeOutIdx t];
                    disp(['taking out trial index ' num2str(t)]);
                    continue;
                end
                wind_max=round(startTimes(t)+window_lift(2));
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
                z_norm_lift_ctx(i,:,:) = (squeeze(r_lift_ctx(i,:,1:2:end))-mu_ctx(i))./(c_soft_z + std_ctx(i));
            end


            %% Decoder
            % frames2predictOver_kine = [15:500];
            % frames2predictOver_neural = [1:500];
            % get data in correct form (movement and neural)
            tmp_kine = traj(selectIdx,:,:); % trials x coord x time
            %kineDecoder = tmp_kine(:,:,frames2predictOver_kine);
            tmp_neural = z_norm_lift_ctx; % neurons x something x time
            tmp_neural = permute(tmp_neural, [2 1 3]);
            save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_DecoderTestMat.mat'],'tmp_kine','tmp_neural');

            %neuralDecoder = tmp_neural(:,:,frames2predictOver_neural);
            % traj indx should be movIdx.nbase
            % define start of lift
%             kineDecoder = tmp_kine(:,:,startFrames_handLift-50:startFrames_handLift+150);
%             neuralDecoder = tmp_neural(:,:,startFrames_handLift-64:startFrames_handLift+150);
            %with APT
            kineDecoder = tmp_kine(:,:,startFrames_APT-50:startFrames_APT+150);
            neuralDecoder = tmp_neural(:,:,startFrames_APT-64:startFrames_APT+150);
            
            %%% kine = 160 x 3 x 2500
            %%% neural = 24 x 130 x 2501
            % Get lift times
            DECODER_FILE = ([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_Decoder.mat']);
            % if ~exist(DECODER_FILE) == 1
            [pred,VAF]=decoder(kineDecoder,neuralDecoder);
            save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_Decoder.mat'],'pred','VAF');
            % else
            %    load(DECODER_FILE);
            %end
            % plot first 9 trials in the z direction
            figure;
            hold on
            for i=1:9
                subplot(3,3,i)
                hold on
                %plot the actual data
                tmp=squeeze(kineDecoder(i,3,:));
                plot(tmp,'r')

                %plot the predicted/decoded value
                tmp=squeeze(pred(i,3,:));
                plot(tmp,'b')
                legend('Observed','Predicted')
                xlabel('Time (ms)')
                ylabel('Position in Z Direction')
            end
            sgtitle([SUB ' ' EXPER_COND ': Variance Explained = ' num2str(VAF)]);
        end
    end
end