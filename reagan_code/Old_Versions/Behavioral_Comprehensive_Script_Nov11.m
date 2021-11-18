% Behavior script - all behavior
clear ALL;
close ALL;
clc;
%% List animal for analysis
animal_idx = 1;
analyze_harmaline = 0; %0=control, 1=harmaline

%% Specify what to plot
plotSampleTrajectories = 0;
%% Add code paths
BASEPATH = 'C:/Users/rcbul/OneDrive - University of North Carolina at Chapel Hill/Hantman_Lab/Harmaline_Project/';
CODE_REAGAN = [BASEPATH 'Code/reagan_code/'];
CODE_CALIB = [BASEPATH 'Code/britton_code/calibration_files/camera_calibration_Jay_7-28-16/Calib_Results_stereo.mat'];
CODE_TRKR = [BASEPATH 'Code/britton_code/tracker'];
CODE_BRITTON = [BASEPATH 'Code/britton_code/code'];
CODE_BRITTON_PLOT = [BASEPATH 'Code/britton_code/other_code'];
addpath(genpath(CODE_REAGAN));
addpath(CODE_CALIB);
addpath(CODE_TRKR);
addpath(genpath(CODE_BRITTON));
addpath(genpath(CODE_BRITTON_PLOT));

Directory_Animals;
%% Prep which experimental condition to analyze
if analyze_harmaline == 1
    exper_condition = 'harm';
else
    exper_condition = 'control';
end
%%
for isub = 1:length(animals)
    % Find how many sessions for animal in experimental condition of choice
    % (control or harmaline)
    SUB = animals{isub};
    ExperSessions = eval(sprintf('%s_%sBehaviorVideos',SUB,exper_condition));
    % Add path for figures
    BEHAV_FIG = [BASEPATH 'Figures\' SUB '\Behavior\'];
    for isession = 1:length(ExperSessions)
        EXPER_SESSION = ExperSessions{isession};
        if isempty(EXPER_SESSION)
            continue;
        end
        TRK = [BASEPATH 'Data/' SUB 'necab1_Chr2/' EXPER_SESSION 'movies/trk'];
        addpath(TRK);
        %% Load Data
        nframes=2500; %variable! %number of frames to analyze. Going to depend on how many frames in movie
        [traj conf] = get_traj_3D_v2(TRK,CODE_CALIB,nframes);
        % traj (trials x coordinate xyz x frames); % x is
        % dim1,vertical
        % conf is the confidence value assigned by the classifier
        
        %% Set trials
        % Define Trials by baseline, stim, washout M340_20210811_controlTrials.nbase
        trialIdxs = eval(sprintf('%s_%s_%sTrials',SUB,EXPER_SESSION,exper_condition));
        base.trialIdxs = trialIdxs.nbase;
        pert.trialIdxs = trialIdxs.npert;
        wash.trialIdxs = trialIdxs.nwash;
        %% Assign trials based on trial type
        % traj(trials x coordinate xyz x frames);
        base.traj = traj(base.trialIdxs,:,:);
        pert.traj = traj(pert.trialIdxs,:,:);
        wash.traj = traj(wash.trialIdxs,:,:);
        %% Set Graph Defaults
        SetGraphDefaults;
        if plotSampleTrajectories == 1
            %% Plot Sample 1D Trajectories in Z
            % INPUT PARSERS AVAILABLE: coordinate_dim, startFrame, endFrame
            plot1DTrajectories(base.traj,'endFrame',nframes);
            sgtitle([SUB ' ' exper_condition ': handpaths baseline trials ZDim']);
            savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' exper_condition '_zDim_baseline.fig']);
            plot1DTrajectories(pert.traj, 'endFrame',nframes);
            sgtitle([SUB ' ' exper_condition ': handpaths stimulation trials ZDim']);
            savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' exper_condition '_zDim_stim.fig']);
            plot1DTrajectories(wash.traj, 'endFrame',nframes);
            sgtitle([SUB ' ' exper_condition ': handpaths washout trials ZDim']);
            savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' exper_condition '_zDim_washout.fig']);
            %% Plot Sample 1D Trajectories in y
            % INPUT PARSERS AVAILABLE: coordinate_dim, startFrame, endFrame
            plot1DTrajectories(base.traj,'endFrame',nframes,'coordinate_dim',2);
            sgtitle([SUB ' ' exper_condition ': handpaths baseline trials YDim']);
            savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' exper_condition '_yDim_baseline.fig']);
            plot1DTrajectories(pert.traj, 'endFrame',nframes,'coordinate_dim',2);
            sgtitle([SUB ' ' exper_condition ': handpaths stimulation trials YDim']);
            savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' exper_condition '_yDim_stim.fig']);
            plot1DTrajectories(wash.traj, 'endFrame',nframes,'coordinate_dim',2);
            sgtitle([SUB ' ' exper_condition ': handpaths washout trials YDim']);
            savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' exper_condition '_yDim_washout.fig']);
            
            %% Plot 3D Trajectories
            plot3DTrajectories(base.traj, 'endFrame',nframes);
            sgtitle([SUB ' ' exper_condition ': 3D handpaths baseline trials']);
            savefig([BASEPATH 'Figures\' SUB '\Behavior\3DTrajectories\' SUB '_' EXPER_SESSION '_' exper_condition '_baseline.fig']);
            plot3DTrajectories(pert.traj, 'endFrame',nframes);
            sgtitle([SUB ' ' exper_condition ': 3D handpaths stimulation trials']);
            savefig([BASEPATH 'Figures\' SUB '\Behavior\3DTrajectories\' SUB '_' EXPER_SESSION '_' exper_condition '_stim.fig']);
            plot3DTrajectories(wash.traj, 'endFrame',nframes);
            sgtitle([SUB ' ' exper_condition ': 3D handpaths washout trials']);
            savefig([BASEPATH 'Figures\' SUB '\Behavior\3DTrajectories\' SUB '_' EXPER_SESSION '_' exper_condition '_washout.fig']);
        end
        %% Find start frame of movement
        % Only for trials with a score of 1 during baseline and stimulation
        [movStrt1.nbase] = getMovStartFrame(1,traj,trialIdxs.trialScore,trialIdxs.nbase);
        [movStrt1.npert] = getMovStartFrame(1,traj,trialIdxs.trialScore,trialIdxs.npert);
        [movStrt1.nwash] = getMovStartFrame(1,traj,trialIdxs.trialScore,trialIdxs.nwash);
        %% Find end frame of reach movement
        % only for trials with a score of 1 suring baseline and stimulation
        [movEnd1.nbase,mov1Idx.nbase] = getMovEndFrame(1,traj,trialIdxs.trialScore,trialIdxs.nbase, movStrt1.nbase);
        [movEnd1.npert,mov1Idx.npert] = getMovEndFrame(1,traj,trialIdxs.trialScore,trialIdxs.npert, movStrt1.npert);
        [movEnd1.nwash,mov1Idx.nwash] = getMovEndFrame(1,traj,trialIdxs.trialScore,trialIdxs.nwash, movStrt1.nwash);
        %% Plot Overlay 1D - only reach
        dimPlot = 3; % if change add dimPlot input to plotOverlay1D function
        % Plot overlay for baseline and stim trials
        plotOverlay1D(movStrt1.nbase, movEnd1.nbase,mov1Idx.nbase,traj);
        title([SUB ' ' exper_condition ': Overlay successful baseline Trials (Dim ' num2str(dimPlot) ')']);
        savefig([BASEPATH 'Figures\' SUB '\Behavior\1DTrajectories\' SUB '_' EXPER_SESSION '_' exper_condition '_overlayReachSuccessBaseline_dim' num2str(dimPlot) '.fig']);
        
        plotOverlay1D(movStrt1.npert, movEnd1.npert,mov1Idx.npert,traj);
        title([SUB ' ' exper_condition ': Overlay successful stim Trials (Dim ' num2str(dimPlot) ')']);
        savefig([BASEPATH 'Figures\' SUB '\Behavior\1DTrajectories\' SUB '_' EXPER_SESSION '_' exper_condition '_overlayReachSuccessStim_dim' num2str(dimPlot) '.fig']);
        
        % PRE and POST baseline
        movStrt1BL = [movStrt1.nbase movStrt1.nwash];
        movEnd1BL = [movEnd1.nbase movEnd1.nwash];
        movStrt1Idxs = [mov1Idx.nbase' mov1Idx.nwash']';
        plotOverlay1D(movStrt1BL, movEnd1BL,movStrt1Idxs,traj);
        title([SUB ' ' exper_condition ': Overlay successful pre and post stim (Dim ' num2str(dimPlot) ')']);
        savefig([BASEPATH 'Figures\' SUB '\Behavior\1DTrajectories\' SUB '_' EXPER_SESSION '_' exper_condition '_overlayReachSuccessAllBL_dim' num2str(dimPlot) '.fig']);
        %% Plot Overlay 1D - whole movement
        % Plot baseline pre stim trials
        plotOverlay1D(movStrt1.nbase, movEnd1.nbase,mov1Idx.nbase, traj,'PlotEntireReach',1);
        title([SUB ' ' exper_condition ': Overlay successful baseline Trials (Dim ' num2str(dimPlot) ')']);
        savefig([BASEPATH 'Figures\' SUB '\Behavior\1DTrajectories\' SUB '_' EXPER_SESSION '_' exper_condition '_overlaySuccessBaseline_dim' num2str(dimPlot) '.fig']);
        % Plot stim
        plotOverlay1D(movStrt1.npert, movEnd1.npert,mov1Idx.npert,traj,'PlotEntireReach',1);
        title([SUB ' ' exper_condition ': Overlay successful stim Trials (Dim ' num2str(dimPlot) ')']);
        savefig([BASEPATH 'Figures\' SUB '\Behavior\1DTrajectories\' SUB '_' EXPER_SESSION '_' exper_condition '_overlaySuccessStim_dim' num2str(dimPlot) '.fig']);
        % Plot baseline pre stim and post stim
        movStrt1BL = [movStrt1.nbase movStrt1.nwash];
        movEnd1BL = [movEnd1.nbase movEnd1.nwash];
        movStrt1Idxs = [mov1Idx.nbase' mov1Idx.nwash']';
        plotOverlay1D(movStrt1BL, movEnd1BL,movStrt1Idxs,traj,'PlotEntireReach',1);
        title([SUB ' ' exper_condition ': Overlay successful pre and post stim (Dim ' num2str(dimPlot) ')']);
        savefig([BASEPATH 'Figures\' SUB '\Behavior\1DTrajectories\' SUB '_' EXPER_SESSION '_' exper_condition '_overlaySuccessAllBL_dim' num2str(dimPlot) '.fig']);
       
        %% Plot Overlay 3D - reach movement
        plotOverlay3D(movStrt1.nbase, movEnd1.nbase,mov1Idx.nbase, traj);
        title([SUB ' ' exper_condition ': Overlay successful baseline Trials']);
        savefig([BASEPATH 'Figures\' SUB '\Behavior\3DTrajectories\' SUB '_' EXPER_SESSION '_' exper_condition '_overlayReachSuccessBaseline.fig']);
        save([BASEPATH 'Data_Analyzed\' SUB '\Behavior\' SUB '_' EXPER_SESSION '_' exper_condition '_Overlay3DVariables.mat'],'movStrt1', 'movEnd1','mov1Idx', 'traj');
        % Plot post and pre stim baseline
        movStrt1BL = [movStrt1.nbase movStrt1.nwash];
        movEnd1BL = [movEnd1.nbase movEnd1.nwash];
        movStrt1Idxs = [mov1Idx.nbase' mov1Idx.nwash']';
        plotOverlay3D(movStrt1BL, movEnd1BL,movStrt1Idxs,traj);
        title([SUB ' ' exper_condition ': Overlay successful pre and post stim']);
        savefig([BASEPATH 'Figures\' SUB '\Behavior\3DTrajectories\' SUB '_' EXPER_SESSION '_' exper_condition '_overlayReachSuccessAllBL.fig']);
       
    end %session
end % animal